#include <NumCpp.hpp>


int main()
{
  nc::NdArray<double> A = nc::random::randFloat<double>(nc::Shape({5, 5}), 0.0, 1.0);
  nc::NdArray<double> B = nc::copy(A);
  // print out for comparison
  std::cout << A << std::endl;
  std::cout << B << std::endl;
  // try column assignment to row
  nc::NdArray<double> C = nc::zeros<double>(5, 5);
  A(0, A.cSlice()) = C(0, C.rSlice());
  // try to use put function on B instead
  B.put(0, B.cSlice(), C(0, C.rSlice()));
  // print out for comparison
  std::cout << A << std::endl;
  std::cout << B << std::endl;

  // for sanity, print type of A, and A(0, 0) and A(0, A.cSlice())
  std::cout << "Type of A: " << typeid(A).name() << std::endl;
  std::cout << "A(0, 0): " << typeid(A(0, 0)).name() << std::endl;
  std::cout << "A(0, A.cSlice()): " << typeid(A(0, A.cSlice())).name() << std::endl;
  return 0;
}
