#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <NumCpp.hpp>

#include <vector>
#include <iterator>
#include <array>
#include <numeric>

// export

#include "splineconcepts_export.h"
#include "Spline.hpp"

class Geometry
{
public:
  nc::NdArray<double> v_;
  nc::NdArray<int> f_;
  nc::NdArray<double> n_;
  nc::NdArray<double> c_;
  nc::NdArray<double> uv_;
  nc::NdArray<double> t_;
  nc::NdArray<double> b_;

public:
  Geometry() = default;
  ~Geometry() = default;

  void fillBetween(const nc::NdArray<double>& pointsLower, const nc:ndArray<double>& pointsUpper)
  {
    if (pointsLower.numCols() != 3 || pointsUpper.numCols() != 3)
    {
      throw std::invalid_argument("Geometry: pointsLower and pointsUpper must have 3 columns.");
    }
    
    const auto A& = (pointsLower.numRows() < pointsUpper.numRows()) ? pointsLower : pointsUpper;
    const auto B& = (pointsLower.numRows() <= pointsUpper.numRows()) ? pointsUpper : pointsLower;
    const auto numA = A.numRows();
    const auto numB = B.numRows();
    const auto numPoints = numA + numB;
    const auto offset = v_.numRows();
    v_.resize(offset + numPoints, 3);
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
        f_.append({offset + Ai, offset + numA + Bi, offset + numA + Bi + 1});
        ++Bi;
      }
      // no matter if an expansion was made, still triangulate the normal pair
      // triangle with two from A and one from B
      f_.append({offset + Ai + 1, offset + Ai, offset + numA + Bi});
      // triangle with one from A and two from B
      f_.append({offset + Ai + 1, offset + numA + Bi, offset + numA + Bi + 1});
      ++Ai;
      ++Bi;
    }
  }

  nc::NdArray<double> kMeans(const nc::NdArray<double>& points, int numClusters = 3)
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
      clusterCenters.put(i, clusterCenters.cSlice(), points(nc::random<int>(0, numPoints - 1), nc::Slice(0, 3)));
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
        newClusterCenters.put(i, newClusterCenters.cSlice(), newClusterCenters(i, nc::Slice(0, 3)) / clusterSizes[i]);
      }
      // check for convergence
      converged = nc::all(nc::all(newClusterCenters == clusterCenters));
      clusterCenters = newClusterCenters;
    }
    return clusterCenters;
  }

  void fillWithin(const nc::NdArray<double>& pointsSurround)
  {
    // persistence based determination how many mid points we need
    const auto numPoints = pointsSurround.numRows();
    const auto offset = v_.numRows();
    double smallest_distance = std::numeric_limits<double>::max();
    for (int i = 0; i < numPoints; ++i)
    {
      const auto distance = nc::norm(pointsSurround(i, nc::Slice(0, 3)) - pointsSurround((i + 1) % numPoints, nc::Slice(0, 3)))[0];
      if (distance < smallest_distance)
      {
        smallest_distance = distance;
      }
    }
    // find the minimal set of points that has direct line of sight to all points
    int numMidpoints = 1;
    std::vector<nc::Vec3> midpoints;
    midpoints.push_back(nc::average(pointsSurround,nc::Axis::ROW));

  }

};


#endif // GEOMETRY_HPP
