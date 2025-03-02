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



/**
 * @brief An interface that defines the methods that a spline class should implement.
 */

class SPLINECONCEPTS_EXPORT Spline
{
  // Construction / Destruction
public:
  virtual ~Spline() = default;

  // Public Methods
public:
  virtual nc::Vec3 operator()(double t) const = 0;
  virtual nc::Vec3 derivative(double t) const = 0;
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
    PREINITIALIZE_SPLINE_VARIABLES();
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

  nc::Vec3 secondDerivative(double t) const
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
    const auto it = std::upper_bound(t_.begin(), t_.end(), t);
    return std::distance(t_.begin(), it);
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
    y_(i, y_.cSlice()) = {point.x, point.y, point.z};
    computeT();
  }

protected:
  void computeT()
  {
    // generate differences
    t_ = nc::diff(y_, nc::Axis::ROW);
    // compute the norm of the differences
    t_ = nc::sqrt(nc::sum(nc::power(t_, 2.0), nc::Axis::COL));
    t_ = nc::cumsum(t_);
    // normalize the differences
    //t_ = t_ / t_.back();

    // compute previous point as extension of the first segment
    prevPoint_ = y_(0, y_.cSlice()) - t_[0] * (y_(1, y_.cSlice()) - y_(0, y_.cSlice()));
    // compute next point as extension of the last segment
    nextPoint_ = y_(-1, y_.cSlice()) + t_[-1] * (y_(-1, y_.cSlice()) - y_(-2, y_.cSlice()));
  }


private:
  nc::NdArray<double> y_;
  nc::NdArray<double> t_;
  nc::Vec3 prevPoint_;
  nc::Vec3 nextPoint_;

  double alpha_;
};



std::shared_ptr<Spline> interpolateSplines(
  std::shared_ptr<Spline> A,
  std::shared_ptr<Spline> B,
  double alpha = 0.5,
  int numPoints = 100)
{
  const auto a = A->toNdArray(numPoints);
  const auto b = B->toNdArray(numPoints);
  const auto c = alpha * a + (1.0 - alpha) * b;
  return std::make_shared<CatmulRomSpline>(c, alpha);
}



#endif // SPLINE_HPP
