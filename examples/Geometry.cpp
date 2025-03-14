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

  return 0;
}