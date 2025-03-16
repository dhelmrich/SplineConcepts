#include <NumCpp.hpp>

#include <Spline.hpp>
#include <RSML.hpp>
#include <Mesh.hpp>
#include <Geometry.hpp>

int main()
{
  // create a geometry object
  Geometry geometry;

  // create some points
  nc::NdArray<double> pointsLower = nc::random::randN<double>({10, 3});
  nc::NdArray<double> pointsUpper = nc::random::randN<double>({10, 3});

  // fill the geometry object
  geometry.fillBetween(pointsLower, pointsUpper);

  // print out the vertices
  std::cout << geometry.v_ << std::endl;

  // print out the triangles
  for (const auto& triangle : geometry.f_)
  {
    std::cout << triangle[0] << " " << triangle[1] << " " << triangle[2] << std::endl;
  }

  return 0;
}