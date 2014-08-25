#include <iomanip>
#include <iostream>
#include "arm.h"

using namespace std;
using namespace gcop;

int main(int argc, char** argv)
{

  Arm arm;

  double p[3];
  
  double a[3] = {0.4, .3, .5};

  cout << "Angles: " << a[0] << ", " << a[1] << ", " << a[2] << endl;

  //  arm.x1 =0;

  arm.Fk(p, a);

  cout << "FK: " << p[0] << ", " << p[1] << ", " << p[2] << endl;

  double as[2][3];
  arm.Ik(as, p);

  cout << "IK0: " << as[0][0] << ", " << as[0][1] << ", " << as[0][2] << endl;
  cout << "IK0: " << as[1][0] << ", " << as[1][1] << ", " << as[1][2] << endl;

  return 0;
}
