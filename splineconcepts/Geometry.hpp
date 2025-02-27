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
    
    auto A = (pointsLower.numRows() < pointsUpper.numRows()) ? pointsLower : pointsUpper;
    auto B = (pointsLower.numRows() <= pointsUpper.numRows()) ? pointsUpper : pointsLower;
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
    v_.slice(nc::Slice(offset, offset + numA)) = A;
    v_.slice(nc::Slice(offset + numA, offset + numPoints)) = B;
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

};


#endif // GEOMETRY_HPP
