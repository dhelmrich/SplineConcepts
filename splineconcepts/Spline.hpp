#ifndef SPLINE_HPP
#define SPLINE_HPP

#define NUMCPP_NO_USE_BOOST

#include <NumCpp.hpp>

#include <vector>
#include <iterator>
#include <array>
#include <numeric>
#include <functional>

// export
#include "splineconcepts_export.h"

// amend hash method for nc::Vec3
namespace std
{
  template <>
  struct hash<nc::Vec3>
  {
    std::size_t operator()(const nc::Vec3& v) const
    {
      return std::hash<double>{}(v.x) ^ std::hash<double>{}(v.y) ^ std::hash<double>{}(v.z);
    }
  };
} // namespace std

class BezierSpline;
class Polynomial;
class Spline;
class CatmulRomSpline;


/**
 * @brief An interface that defines the methods that a spline class should implement.
 */

class SPLINECONCEPTS_EXPORT Spline
{
  // Construction / Destruction
public:
  virtual ~Spline() = default;
  bool intersects(const Spline& other) const
  {
    // coarse resolution of finding sampling over t \in [0, 1]^2
    const int numPoints = 10;
    for (int i = 0; i < numPoints; ++i)
    {
      const double t0 = static_cast<double>(i) / (numPoints - 1);
      for (int j = 0; j < numPoints; ++j)
      {
        const double t1 = static_cast<double>(j) / (numPoints - 1);
        if ((operator()(t0) - other(t1)).norm() < 1e-6)
        {
          return true;
        }
      }
    }
    return false;
  }
  std::vector<std::pair<double,double>> intersection(const Spline& other) const
  {
    std::vector<std::pair<double,double>> result;
    // coarse resolution of finding sampling over t \in [0, 1]^2
    const int numPoints = 10;
    for (int i = 0; i < numPoints; ++i)
    {
      const double t0 = static_cast<double>(i) / (numPoints - 1);
      for (int j = 0; j < numPoints; ++j)
      {
        const double t1 = static_cast<double>(j) / (numPoints - 1);
        if ((operator()(t0) - other(t1)).norm() < 1e-6)
        {
          result.push_back(std::make_pair(t0, t1));
        }
      }
    }
    return result;
  }
  // Public Methods
public:
  virtual nc::Vec3 operator()(double t) const = 0;
  virtual nc::Vec3 derivative(double t) const = 0;
  virtual nc::Vec3 secondDerivative(double t) const = 0;
  virtual nc::NdArray<double> toNdArray(const int numPoints) const = 0;
  virtual double length(double t0, double t1) const = 0;
  virtual double length() const = 0;
  virtual nc::rotations::Quaternion orientation(double t) const = 0;
  virtual std::size_t hash() const = 0;
};

class SPLINECONCEPTS_EXPORT CatmulRomSpline : public Spline
{
public:
  // this class always contains technically more splines than one, but there is never a use for a single spline
  // So we will build the math as long matrices per section
  CatmulRomSpline(nc::NdArray<double> controlPoints, double alpha = 0.5)
      : y_(controlPoints), alpha_(alpha)
  {
    if (y_.numRows() < 4)
    {
      throw std::invalid_argument("CatmulRomSpline: controlPoints must have at least 4 rows.");
    }
    if (y_.numCols() != 3)
    {
      throw std::invalid_argument("CatmulRomSpline: controlPoints must have 3 columns.");
    }
    computeT();
  }

  virtual ~CatmulRomSpline() = default;

#define PREINITIALIZE_SPLINE_VARIABLES() \
  const auto segment = findSegment(t); \
  const auto tSegment = (t - t_[segment]) / (t_[segment + 1] - t_[segment]); \
  const nc::Vec3 p0 = (segment == 0) ? prevPoint_ : y_(segment - 1, y_.cSlice()); \
  const nc::Vec3 p1 = y_(segment, y_.cSlice()); \
  const nc::Vec3 p2 = y_(segment + 1, y_.cSlice()); \
  const nc::Vec3 p3 = (segment == y_.numRows() - 2) ? nextPoint_ : y_(segment + 2, y_.cSlice());

  nc::Vec3 operator()(double t) const override
  {
    const auto segment = findSegment(t);
    const auto tSegment = (t - t_[segment]) / (t_[segment + 1] - t_[segment]);
    const nc::Vec3 p0 = (segment == 0) ? prevPoint_ : y_(segment - 1, y_.cSlice());
    const nc::Vec3 p1 = y_(segment, y_.cSlice());
    const nc::Vec3 p2 = y_(segment + 1, y_.cSlice());
    const nc::Vec3 p3 = (segment == y_.numRows() - 2) ? nextPoint_ : y_(segment + 2, y_.cSlice());
    const nc::Vec3 quadp = 0.5 * ((2.0 * p1) + (-p0 + p2) * tSegment + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * tSegment * tSegment + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * tSegment * tSegment * tSegment);
    const nc::Vec3 linp = p1 + tSegment * (p2 - p1);
    return quadp * alpha_ + linp * (1.0 - alpha_);
  }

  nc::Vec3 derivative(double t) const override
  {
    PREINITIALIZE_SPLINE_VARIABLES();
    const nc::Vec3 quadp = 0.5 * ((-p0 + p2) + 2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * tSegment + 3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * tSegment * tSegment);
    const nc::Vec3 linp = p2 - p1;
    return quadp * alpha_ + linp * (1.0 - alpha_);
  }

  nc::Vec3 secondDerivative(double t) const override
  {
    PREINITIALIZE_SPLINE_VARIABLES();
    const nc::Vec3 quadp = 0.5 * (2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) + 6.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * tSegment);
    return quadp * alpha_;
  }

  nc::NdArray<double> toNdArray(const int numPoints) const override
  {
    nc::NdArray<double> result = nc::zeros<double>(numPoints, 3);
    const double dt = 1.0 / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i)
    {
      result.put(i, result.cSlice(), operator()(i * dt).toNdArray());
    }
    return result;
  }

  std::size_t hash() const override
  {
    std::size_t hash = 0;
    for (int i = 0; i < y_.numRows(); ++i)
    {
      hash += std::hash<nc::Vec3>{}(y_(i, y_.cSlice()));
    }
    return hash;
  }

  double length(double t0, double t1) const override
  {
    if (t0 > t1)
    {
      throw std::invalid_argument("CatmulRomSpline: t0 must be less than or equal to t1.");
    }
    if (t0 < 0.0 || t1 > t_.back())
    {
      throw std::invalid_argument("CatmulRomSpline: t0 and t1 must be in the range [0, 1].");
    }

    const auto segment0 = findSegment(t0);
    const auto segment1 = findSegment(t1);

    double length = 0.0;
    for (int i = segment0; i < segment1; ++i)
    {
      const double tStart = (i == segment0) ? t0 : t_[i];
      const double tEnd = (i == segment1) ? t1 : t_[i + 1];
      const auto p0 = (i == 0) ? prevPoint_ : y_(i - 1, y_.cSlice());
      const auto p1 = y_(i, y_.cSlice());
      const auto p2 = y_(i + 1, y_.cSlice());
      const auto p3 = (i == y_.numRows() - 2) ? nextPoint_ : y_(i + 2, y_.cSlice());

      const auto a = 6.0 * p0 - 12.0 * p1 + 6.0 * p2;
      const auto b = -3.0 * p0 + 9.0 * p1 - 9.0 * p2 + 3.0 * p3;
      const auto c = 3.0 * p1 - 3.0 * p0;

      const auto tStart2 = tStart * tStart;
      const auto tStart3 = tStart2 * tStart;
      const auto tEnd2 = tEnd * tEnd;
      const auto tEnd3 = tEnd2 * tEnd;

      length += (1.0 / 6.0) * (a * (tEnd3 - tStart3) + b * (tEnd2 - tStart2) + c * (tEnd - tStart)).norm();
    }

    return length;
  }

  double length() const override
  {
    return length(0.0, 1.0);
  }

  nc::rotations::Quaternion orientation(double t) const override
  {
    const nc::Vec3 d = derivative(t);
    const nc::Vec3 dd = secondDerivative(t);
    const nc::Vec3 up = nc::Vec3::up();
    const nc::Vec3 right = nc::Vec3::right();
    const nc::Vec3 forward = nc::Vec3::forward();
    // use the derivative matrix and the global coordinate system to compute the orientation
    nc::Vec3 normal = dd.cross(d).normalize();
    const nc::Vec3 binormal = d.cross(normal).normalize();
    const nc::Vec3 tangent = d.normalize();
    const auto dcm = nc::NdArray<double>({right.dot(tangent), right.dot(binormal), right.dot(normal),
                                          up.dot(tangent), up.dot(binormal), up.dot(normal),
                                          forward.dot(tangent), forward.dot(binormal), forward.dot(normal)})
                         .reshape(3, 3);
    return nc::rotations::Quaternion(dcm);
  }

  inline int findSegment(double t) const
  {
    double t0 = 0.0;
    for (int i = 0; i < t_.size(); ++i)
    {
      if (t_[i] > t)
      {
        return i - 1;
      }
    }
    return t_.size() - 2;
  }

  nc::Vec3 getControlPoint(int i) const
  {
    if(i == -1) return prevPoint_;
    if(i == y_.numRows()) return nextPoint_;
    if(i < -1 || i >= y_.numRows()) throw std::invalid_argument("CatmulRomSpline: index out of bounds.");
    return y_(i, y_.cSlice());
  }

  void setControlPoint(int i, const nc::Vec3& point)
  {
    y_.put(i, y_.cSlice(),{point.x, point.y, point.z});
    computeT();
  }

protected:
  void computeT()
  {
    const auto n = y_.numRows();
    // generate differences
    t_ = nc::diff(y_, nc::Axis::ROW);
    // compute the norm of the differences
    t_ = nc::sqrt(nc::sum(nc::power(t_, 2.0), nc::Axis::COL));
    t_ = nc::cumsum(t_);
    // add zero at the beginning
    t_ = nc::insert(t_, 0, 0.0, nc::Axis::COL);

    // compute previous point as extension of the first segment
    const nc::Vec3 diff1 = y_(1, y_.cSlice()) - y_(0, y_.cSlice());
    prevPoint_ = y_(0, y_.cSlice());
    prevPoint_ = prevPoint_ - diff1;
    // compute next point as extension of the last segment
    const nc::Vec3 diff2 = y_(n-1, y_.cSlice()) - y_(n-2, y_.cSlice());
    nextPoint_ = y_(n-1, y_.cSlice());
    nextPoint_ = nextPoint_ + diff2;
  }


private:
  nc::NdArray<double> y_;
  nc::NdArray<double> t_;
  nc::Vec3 prevPoint_;
  nc::Vec3 nextPoint_;

  double alpha_;
};

class SPLINECONCEPTS_EXPORT BezierSpline : public Spline
{
  friend class Polynomial;
public:
  // constructor using control points
  BezierSpline(const nc::NdArray<double>& points, const nc::NdArray<double>& controlPoints)
  {
    if (points.numRows() != controlPoints.numRows() + 1)
    {
      throw std::invalid_argument("BezierSpline: controlPoints must have one less row than points.");
    }
    if (points.numCols() != 3 || controlPoints.numCols() != 3)
    {
      throw std::invalid_argument("BezierSpline: points and controlPoints must have 3 columns.");
    }
    if (points.numRows() < 2)
    {
      throw std::invalid_argument("BezierSpline: points must have at least 2 rows.");
    }
    y_ = points;
    c_ = controlPoints;
    computeT();
  }

  // constructor using tangents
  BezierSpline(const nc::NdArray<double>& points, const nc::NdArray<double>& tangents, bool useTangents = true)
  {
    if (points.numRows() != tangents.numRows())
    {
      throw std::invalid_argument("BezierSpline: points and tangents must have the same number of rows.");
    }
    if (points.numCols() != 3 || tangents.numCols() != 3)
    {
      throw std::invalid_argument("BezierSpline: points and tangents must have 3 columns.");
    }
    if (points.numRows() < 2)
    {
      throw std::invalid_argument("BezierSpline: points must have at least 2 rows.");
    }
    y_ = points;
    // compute control points from tangents
    c_ = nc::zeros<double>(y_.numRows(), 3);
    for (int i = 0; i < y_.numRows(); ++i)
    {
      if (i == 0)
      {
        c_.put(i, c_.cSlice(), y_(i, y_.cSlice()) + tangents(i, tangents.cSlice()) / 3.0);
      }
      else if (i == y_.numRows() - 1)
      {
        c_.put(i, c_.cSlice(), y_(i, y_.cSlice()) - tangents(i - 1, tangents.cSlice()) / 3.0);
      }
      else
      {
        c_.put(i, c_.cSlice(), y_(i, y_.cSlice()) - tangents(i - 1, tangents.cSlice()) / 3.0);
      }
    }
  }

  virtual ~BezierSpline() = default;


public:

#define PREINITIALIZE_BEZIER_SPLINE_VARIABLES() \
  const auto segment = findSegment(t); \
  const auto tSegment = (t - t_[segment]) / (t_[segment + 1] - t_[segment]); \
  const nc::Vec3 p0 = y_(segment    , y_.cSlice()); \
  const nc::Vec3  c = c_(segment    , c_.cSlice()); \
  const nc::Vec3 p1 = y_(segment + 1, y_.cSlice());

  nc::Vec3 operator()(double t) const override
  {
    PREINITIALIZE_BEZIER_SPLINE_VARIABLES();
    return (1.0 - tSegment) * (1.0 - tSegment) * p0 + 2.0 * (1.0 - tSegment) * tSegment * c + tSegment * tSegment * p1;
  }

  nc::Vec3 normal(double t) const
  {
    const auto d = derivative(t);
    const auto dd = secondDerivative(t);
    return dd.cross(d).normalize();
  }

  nc::Vec3 derivative(double t) const override
  {
    PREINITIALIZE_BEZIER_SPLINE_VARIABLES();
    return 2.0 * (1.0 - tSegment) * (c - p0) + 2.0 * tSegment * (p1 - c);
  }

  nc::Vec3 secondDerivative(double t) const
  {
    PREINITIALIZE_BEZIER_SPLINE_VARIABLES();
    return 2.0 * (p1 - 2.0 * c + p0);
  }

  nc::Vec3 deCasteljau(double t) const
  {
    const auto n = y_.numRows() - 1;
    std::vector<nc::Vec3> points;
    points.reserve(n + 1);
    for (int i = 0; i <= n; ++i)
    {
      points.push_back(y_(i, y_.cSlice()));
    }
    for (int i = 1; i <= n; ++i)
    {
      for (int j = 0; j <= n - i; ++j)
      {
        points[j] = (1.0 - t) * points[j] + t * points[j + 1];
      }
    }
    return points[0];
  }

  nc::NdArray<double> toNdArray(const int numPoints) const override
  {
    nc::NdArray<double> result = nc::zeros<double>(numPoints, 3);
    const double dt = 1.0 / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i)
    {
      result.put(i, result.cSlice(), operator()(i * dt).toNdArray());
    }
    return result;
  }

  double length(double t0, double t1) const override
  {
    // for the segment length, use the linear interpolation property of the Bezier curve
    const auto segment0 = findSegment(t0);
    const auto segment1 = findSegment(t1);
    double length = 0.0;
    for (int segment = segment0; segment < segment1; ++segment)
    {
      const auto p0 = y_(segment, y_.cSlice());
      const auto cc = c_(segment, c_.cSlice());
      const auto p1 = y_(segment + 1, y_.cSlice());
      const auto a = 2.0 * (cc - p0);
      const auto b = 2.0 * (p1 - cc);
      const auto tStart = (segment == segment0) ? t0 : t_[segment];
      const auto tEnd = (segment == segment1) ? t1 : t_[segment + 1];
      const auto tStart2 = tStart * tStart;
      const auto tEnd2 = tEnd * tEnd;
      length += (1.0 / 6.0) * nc::norm((a * (tEnd2 - tStart2) + b * (tEnd2 - tStart2)))[0];
    }
    return length;
  }

  double length() const override
  {
    return length(0.0, 1.0);
  }

  void insertControlPoints(double t)
  {
    const auto segment = findSegment(t);
    const auto tSegment = (t - t_[segment]) / (t_[segment + 1] - t_[segment]);
    const auto p0 = y_(segment, y_.cSlice());
    const auto p1 = y_(segment + 1, y_.cSlice());
    const auto cc = (1.0 - tSegment) * p0 + tSegment * p1;
    y_ = nc::insert(y_, segment + 1, cc, nc::Axis::ROW);
    c_ = nc::insert(c_, segment, cc, nc::Axis::ROW);
    computeT();
  }

  nc::rotations::Quaternion orientation(double t) const override
  {
    const nc::Vec3 d = derivative(t);
    const nc::Vec3 dd = secondDerivative(t);
    const nc::Vec3 up = nc::Vec3::up();
    const nc::Vec3 right = nc::Vec3::right();
    const nc::Vec3 forward = nc::Vec3::forward();
    // use the derivative matrix and the global coordinate system to compute the orientation
    nc::Vec3 normal = dd.cross(d).normalize();
    const nc::Vec3 binormal = d.cross(normal).normalize();
    const nc::Vec3 tangent = d.normalize();
    const auto dcm = nc::NdArray<double>({right.dot(tangent), right.dot(binormal), right.dot(normal),
                                          up.dot(tangent), up.dot(binormal), up.dot(normal),
                                          forward.dot(tangent), forward.dot(binormal), forward.dot(normal)})
                         .reshape(3, 3);
    return nc::rotations::Quaternion(dcm);
  }

  std::size_t hash() const override
  {
    std::size_t hash = 0;
    for (int i = 0; i < y_.numRows(); ++i)
    {
      hash += std::hash<nc::Vec3>{}(y_(i, y_.cSlice()));
    }
    return hash;
  }

  const nc::NdArray<double>& getPoints() const
  {
    return y_;
  }

  nc::NdArray<double> getControlPoints(int AlternativeResolution) const
  {
    nc::NdArray<double> altControlPoints = nc::zeros<double>(AlternativeResolution - 1, 3);
    for (int i = 0; i < AlternativeResolution - 1; ++i)
    {
      const double t = static_cast<double>(i) / (AlternativeResolution - 1);
      const auto segment = findSegment(t);
      const auto tSegment = (t - t_[segment]) / (t_[segment + 1] - t_[segment]);
      const nc::NdArray<double> p0 = y_(segment, y_.cSlice());
      const nc::NdArray<double> cc = c_(segment, c_.cSlice());
      const nc::NdArray<double> p1 = y_(segment + 1, y_.cSlice());
      altControlPoints.put(i, altControlPoints.cSlice(), (1.0 - tSegment) * cc + tSegment * p1);
    }
    return std::move(altControlPoints);
  }

  const nc::NdArray<double>& getControlPoints() const
  {
    return c_;
  }

  const nc::NdArray<double>& getT() const
  {
    return t_;
  }

  


  private:
  nc::NdArray<double> y_;
  nc::NdArray<double> c_; // control points from tangents
  nc::NdArray<double> t_;

  void computeT()
  {
    // the t values are the distances between the points
    t_ = nc::diff(y_, nc::Axis::ROW);
    t_ = nc::sqrt(nc::sum(nc::power(t_, 2.0), nc::Axis::COL));
    t_ = nc::cumsum(t_);
  }

  inline int findSegment(double t) const
  {
    const auto it = std::upper_bound(t_.begin(), t_.end(), t);
    return std::distance(t_.begin(), it);
  }
};



template <typename SplineType>
std::shared_ptr<SplineType> interpolateSplines(
  std::shared_ptr<SplineType> A,
  std::shared_ptr<SplineType> B,
  double alpha = 0.5,
  int numPoints = 100)
{
  static_assert(std::is_base_of<Spline, SplineType>::value, "SplineType must derive from Spline");

  const auto a = A->toNdArray(numPoints);
  const auto b = B->toNdArray(numPoints);
  const auto c = alpha * a + (1.0 - alpha) * b;
  return std::make_shared<SplineType>(c);
}

template <>
std::shared_ptr<BezierSpline> interpolateSplines(
  std::shared_ptr<BezierSpline> A,
  std::shared_ptr<BezierSpline> B,
  double alpha,
  int numPoints)
{
  const auto aPoints = A->getPoints();
  const auto bPoints = B->getPoints();
  const auto aTangents = A->getControlPoints();
  const auto bTangents = B->getControlPoints();
  
  const auto cPoints = alpha * aPoints + (1.0 - alpha) * bPoints;
  const auto cTangents = alpha * aTangents + (1.0 - alpha) * bTangents;

  return std::make_shared<BezierSpline>(cPoints, cTangents, true);
}


/**
 * @brief A class that represents a polynomial in the form of a0 + a1 * x + a2 * x^2 + ... + an * x^n.
 */
class SPLINECONCEPTS_EXPORT Polynomial
{
public:
  // Construction / Destruction
  Polynomial(const nc::NdArray<double>& coefficients)
      : coefficients_(coefficients)
  {
  }

  Polynomial() = default;

  // Public Methods
public:
  nc::NdArray<double> operator()(const nc::NdArray<double>& x) const
  {
    nc::NdArray<double> result = nc::zeros<double>(x.shape());
    for (int i = 0; i < coefficients_.size(); ++i)
    {
      result += coefficients_[i] * nc::power(x, i);
    }
    return result;
  }

  Polynomial derivative() const
  {
    nc::NdArray<double> newCoefficients = nc::zeros<double>(coefficients_.shape());
    for (int i = 1; i < coefficients_.size(); ++i)
    {
      newCoefficients[i - 1] = i * coefficients_[i];
    }
    return Polynomial(newCoefficients);
  }

  Polynomial integral() const
  {
    nc::NdArray<double> newCoefficients = nc::zeros<double>(coefficients_.shape());
    for (int i = 0; i < coefficients_.size(); ++i)
    {
      newCoefficients[i + 1] = coefficients_[i] / (i + 1);
    }
    return Polynomial(newCoefficients);
  }

  inline static int factorial(int n)
  {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
  }

  inline static int nCr(int n, int r)
  {
    if (r == 0 || r == n)
    {
      return 1;
    }
    return nCr(n - 1, r) + nCr(n - 1, r - 1);
  }

  static Polynomial fromBezier(const BezierSpline& bezier)
  {
    Polynomial result;
    const auto y = bezier.y_;
    const auto c = bezier.c_;
    const auto n = y.numRows() - 1;
    nc::NdArray<double> coefficients = nc::zeros<double>(n + 1);
    for (int i = 0; i < n + 1; ++i)
    {
      coefficients[i] = 0.0;
      for (int j = 0; j < i + 1; ++j)
      {
        coefficients[i] += nCr(i, j) * y(i, 0) * nc::power(1.0 - c(i, 0), i - j) * nc::power(c(i, 0), j);
      }
    }
    return result;
  }

private:
  nc::NdArray<double> coefficients_;
};

enum class CouplineType
{
  FIXED, // if the coupling point should not move under any circumstances, even when scaling the splines
  SECTION, // if the coupling point is relative to the spline coordinates
  POINT, // if the coupling point is relative to the global coordinates
  STRETCH, // same as SECTION, but forces the splines to scale together

};

struct SPLINECONCEPTS_EXPORT SplineCoupling
{

};

/**
 * @brief A class that represents a group of splines that are coupled together.
 * This group can describe a 3D or even 2D shape, and is the expression of how splines form a mesh section.
 * The coupling is the enforcing of transformation rules, and the group is the collection of splines.
 * User-facing: Operations that are targeted at a spline in a group should by default instead target the group.
 * Coupling: By default, splines in a group are coupled using STRETCH coupling.
 */
class SPLINECONCEPTS_EXPORT Group
{

};

#endif // SPLINE_HPP
