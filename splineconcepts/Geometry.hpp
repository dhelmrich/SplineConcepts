#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <NumCpp.hpp>

#include <vector>
#include <iterator>
#include <array>
#include <optional>
#include <algorithm>
#include <functional>
#include <tuple>
#include <numeric>

// export

#include "splineconcepts_export.h"
#include "Spline.hpp"

class Geometry
{
public:
  nc::NdArray<double> v_;
  std::vector<std::array<nc::NdArray<double>::size_type, 3>> f_;
  nc::NdArray<double> n_;
  nc::NdArray<double> c_;
  nc::NdArray<double> uv_;
  nc::NdArray<double> t_;
  nc::NdArray<double> b_;

public:
  Geometry() = default;
  ~Geometry() = default;

  void fillBetween(const nc::NdArray<double>& pointsLower, const nc::NdArray<double>& pointsUpper)
  {
    if (pointsLower.numCols() != 3 || pointsUpper.numCols() != 3)
    {
      throw std::invalid_argument("Geometry: pointsLower and pointsUpper must have 3 columns.");
    }
    
    const auto& A = (pointsLower.numRows() < pointsUpper.numRows()) ? pointsLower : pointsUpper;
    const auto& B = (pointsLower.numRows() <= pointsUpper.numRows()) ? pointsUpper : pointsLower;
    const auto numA = A.numRows();
    const auto numB = B.numRows();
    const auto numPoints = numA + numB;
    const auto offset = v_.numRows();
    v_ = nc::zeros<double>(numPoints, 3);
    // if there is a difference in the number of points, we choose evenly spaced indices in the smaller array that expand to the larger array
    int diff = numB - numA;
    auto indices = nc::linspace<int>(0, numA - 1, numPoints);
    int Ai = 0;
    int Bi = 0;
    // prefill the vertices
    v_.put(nc::Slice(offset, offset + numA), v_.cSlice(), A);
    v_.put(nc::Slice(offset + numA, offset + numPoints), v_.cSlice(), B);
    // fill the triangles
    while (Ai < numA && Bi < numB)
    {
      if(indices.contains(Ai))
      {
        // triangle with two from B and one from A
        f_.push_back({offset + Ai, offset + numA + Bi, offset + numA + Bi + 1});
      }
      // no matter if an expansion was made, still triangulate the normal pair
      // triangle with two from A and one from B
      f_.push_back({offset + Ai + 1, offset + numA + Bi, offset + numA + Bi + 1});
      // triangle with one from A and two from B
      f_.push_back({offset + Ai, offset + Ai + 1, offset + numA + Bi + 1});
      ++Ai;
      ++Bi;
    }
  }

  std::tuple<nc::NdArray<double>,nc::NdArray<int>> kMeans(const nc::NdArray<double>& points, int numClusters = 3, double* rating = nullptr)
  {
    if (points.numCols() != 3)
    {
      throw std::invalid_argument("Geometry: points must have 3 columns.");
    }
    const auto numPoints = points.numRows();
    // initialize the cluster centers
    nc::NdArray<double> clusterCenters = nc::zeros<double>(numClusters, 3);
    for (int i = 0; i < numClusters; ++i)
    {
      clusterCenters.put(i, clusterCenters.cSlice(), points(nc::random::randInt<int>(0, numPoints - 1), nc::Slice(0, 3)));
    }
    // initialize the cluster assignments
    nc::NdArray<int> clusterAssignments = nc::zeros<int>(numPoints);
    // iterate until convergence
    bool converged = false;
    while (!converged)
    {
      // assign each point to the nearest cluster
      for (int i = 0; i < numPoints; ++i)
      {
        const auto point = points(i, nc::Slice(0, 3));
        double minDistance = std::numeric_limits<double>::max();
        int minCluster = 0;
        for (int j = 0; j < numClusters; ++j)
        {
          const auto clusterCenter = clusterCenters(j, nc::Slice(0, 3));
          const auto distance = nc::norm(point - clusterCenter)[0];
          if (distance < minDistance)
          {
            minDistance = distance;
            minCluster = j;
          }
        }
        clusterAssignments[i] = minCluster;
      }
      // update the cluster centers
      nc::NdArray<double> newClusterCenters = nc::zeros<double>(numClusters, 3);
      nc::NdArray<int> clusterSizes = nc::zeros<int>(numClusters);
      for (int i = 0; i < numPoints; ++i)
      {
        const auto cluster = clusterAssignments[i];
        newClusterCenters.put(cluster, newClusterCenters.cSlice(), newClusterCenters(cluster, nc::Slice(0, 3)) + points(i, nc::Slice(0, 3)));
        ++clusterSizes[cluster];
      }
      for (int i = 0; i < numClusters; ++i)
      {
        newClusterCenters.put(i, newClusterCenters.cSlice(), newClusterCenters(i, nc::Slice(0, 3)) / static_cast<double>(clusterSizes[i]));
      }
      // check for convergence
      converged = nc::all(nc::all(newClusterCenters == clusterCenters)).item();
      clusterCenters = newClusterCenters;
    }
    // calculate the rating if requested
    if (rating != nullptr)
    {
      *rating = 0.0;
      for (int i = 0; i < numPoints; ++i)
      {
        const auto point = points(i, nc::Slice(0, 3));
        const auto clusterCenter = clusterCenters(clusterAssignments[i], nc::Slice(0, 3));
        *rating += nc::norm(point - clusterCenter)[0];
      }
    }
    return std::make_tuple(clusterCenters, clusterAssignments);
  }

  void fillWithin(const nc::NdArray<double>& pointsSurround, std::optional<int> fixedMidpoints = std::nullopt)
  {
    // persistence based determination how many mid points we need
    const auto numPoints = pointsSurround.numRows();
    const auto offset = v_.numRows();
    double rating = 0.0;
    double ratingPrev = std::numeric_limits<double>::max();
    // find the minimal set of points that has direct line of sight to all points
    int numMidpoints = 1;
    nc::NdArray<double> midpoints;
    nc::NdArray<int> clusterAssignments;
    if (fixedMidpoints.has_value())
    {
      numMidpoints = fixedMidpoints.value();
    }
    else while(rating < ratingPrev)
    {
      ratingPrev = rating;
      tie(midpoints, clusterAssignments) = kMeans(pointsSurround, numMidpoints, &rating);
    }
    // fill the vertices
    v_.reshape(offset + numPoints + numMidpoints, 3);
    v_.put(nc::Slice(offset, offset + numPoints), v_.cSlice(), pointsSurround);
    v_.put(nc::Slice(offset + numPoints, offset + numPoints + numMidpoints), v_.cSlice(), midpoints);
    // fill the triangles
    for (int i = 0; i < numPoints; ++i)
    {
      f_.push_back({offset + i, offset + (i + 1) % numPoints, offset + numPoints + clusterAssignments[i]});
    }
    
  }

  void extrudeAround(const Spline& spline, const nc::NdArray<double>& base, int resolution = 20, std::optional<nc::NdArray<double>> radius = std::nullopt)
  {
    if (base.numCols() != 3)
    {
      throw std::invalid_argument("Geometry: base must have 3 columns.");
    }
    if (radius.has_value() && radius.value().numCols() != 1)
    {
      throw std::invalid_argument("Geometry: radius must have 1 column.");
    }
    // base: n x 3, to be placed perpendicularly around the spline, triangles filling between base shapes
    // resolution: number of times the base shape is placed around the spline
    // radius: n x 1, radius of the base shape
    const auto numBase = base.numRows();
    const auto numVertices = numBase * resolution;
    const auto offset = v_.numRows();
    v_ = nc::zeros<double>(numVertices, 3);
    // fill the geometry
    for (int i = 0; i < resolution; ++i)
    {
      const auto t = static_cast<double>(i) / (resolution - 1);
      const auto point = spline(t);
      const auto tangent = spline.derivative(t);
      const auto normal = spline.secondDerivative(t);
      const auto binormal = tangent.cross(normal).normalize();
      for (int j = 0; j < numBase; ++j)
      {
        const auto radiusValue = (radius.has_value()) ? (radius.value())[j] : 1.0;
        // rotate the base shape with the quaternion
        const auto basePoint = base(j, nc::Slice(0, 3));
        // place the base shape around the spline
        const auto finalPoint = point + radiusValue * (basePoint[0] * normal + basePoint[1] * binormal + basePoint[2] * tangent);
        v_.put(offset + i * numBase + j, v_.cSlice(), finalPoint.toNdArray());
        // append triangles for this point
        if (i > 0 && j > 0)
        {
          f_.push_back({offset + i * numBase + j, offset + i * numBase + j - 1, offset + (i - 1) * numBase + j - 1});
          f_.push_back({offset + i * numBase + j, offset + (i - 1) * numBase + j - 1, offset + (i - 1) * numBase + j});
        }
        // normal = point - spline point
        n_.put(offset + i * numBase + j, n_.cSlice(), (finalPoint - point).normalize().toNdArray());
        // texture coordinates: u = t, v = j / numBase
        uv_.put(offset + i * numBase + j, uv_.cSlice(), {t, static_cast<double>(j) / numBase});
        // tangent = spline derivative
        t_.put(offset + i * numBase + j, t_.cSlice(), tangent.toNdArray());
        
      }
    }
  }

};


#endif // GEOMETRY_HPP
