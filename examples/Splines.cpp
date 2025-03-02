#include <vector>
#include <iterator>
#include <iostream>
#include <array>
#include "Spline.hpp"



// A test that makes sure that the spline can do a straight line

int main()
{
  nc::NdArray<double> c = nc::zeros<double>(10, 3);
  c.put(c.rSlice(), 0, nc::linspace<double>(0, 9, 10, nc::EndPoint::YES));
  const auto spline = std::make_shared<CatmulRomSpline>(c, 0);

  std::cout << "Spline points: " << std::endl;
  for(int i = -1; i < 11; ++i)
  {
    std::cout << spline->getControlPoint(i) << " from " << c(i, c.cSlice()) << std::endl;
  }

  for (int i = 0; i < 10; ++i)
  {
    const double t = static_cast<double>(i) / 9.0;
    const auto point = spline->operator()(t);
    std::cout << "spline(" << t << ") = " << point << std::endl;
  }
}
