#include <vector>
#include <iterator>
#include <array>

#include "Spline.hpp"


#include <gtest/gtest.h>

#ifndef SPLINE_TEST_VERBOSE
#define SPLINE_TEST_VERBOSE 0
#endif


// A test to make sure findSegment works
TEST(Spline, FindSegment)
{
  const nc::NdArray<double> controlPoints = nc::append(nc::linspace<double>(0, 10, 10).reshape(10, 1), nc::zeros<double>(10, 2), nc::Axis::COL);
  
  //const auto spline = std::make_shared<CatmulRomSpline>(controlPoints);
  // if we have a verbose test, print out the control points
  if (SPLINE_TEST_VERBOSE)
  {
    controlPoints.print();

    // flush output
    std::cout << std::endl;
    
    // generate differences
    auto t_ = nc::diff(controlPoints, nc::Axis::ROW);
    std::cout << t_.str() << std::endl;
    // compute the norm of the differences
    t_ = nc::sqrt(nc::sum(nc::power(t_, 2.0), nc::Axis::COL));
    std::cout << t_.str() << std::endl;
    // normalize the differences
    t_ = nc::cumsum(t_);
    std::cout << t_.str() << std::endl;
    t_ = t_ / t_.back();
    std::cout << t_.str() << std::endl;

    // compute previous point as extension of the first segment
    auto prevPoint_ = controlPoints(0, controlPoints.cSlice()) - t_[0] * (controlPoints(1, controlPoints.cSlice()) - controlPoints(0, controlPoints.cSlice()));
    std::cout << prevPoint_.str() << std::endl;

    auto nextPoint_ = controlPoints(-1, controlPoints.cSlice()) + t_[-1] * (controlPoints(-1, controlPoints.cSlice()) - controlPoints(-2, controlPoints.cSlice()));
    std::cout << nextPoint_.str() << std::endl;
  }

  const auto spline = std::make_shared<CatmulRomSpline>(controlPoints, 0);
  EXPECT_DOUBLE_EQ(spline->findSegment(0.0), 0);
  EXPECT_DOUBLE_EQ(spline->findSegment(0.5), 0);
  EXPECT_DOUBLE_EQ(spline->findSegment(1.0), 1);
  
}

// A test that makes sure that the spline can do a straight line

TEST(Spline, StraightLine)
{
  nc::NdArray<double> controlPoints = nc::zeros<double>(10, 3);
  controlPoints(controlPoints.rSlice(), 0) = nc::linspace<double>(0, 10, 10);
  const auto spline = std::make_shared<CatmulRomSpline>(controlPoints, 0);

  for (int i = 0; i < 10; ++i)
  {
    const double t = static_cast<double>(i) / 9.0;
    const auto point = spline->operator()(t);
    EXPECT_DOUBLE_EQ(point.x, t * 10.0);
    EXPECT_DOUBLE_EQ(point.y, 0.0);
    EXPECT_DOUBLE_EQ(point.z, 0.0);
  }
}
