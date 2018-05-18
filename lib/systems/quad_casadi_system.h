#ifndef QUADCASADISYSTEM_H
#define QUADCASADISYSTEM_H
#include "casadi_system.h"
#include "rn.h"

namespace gcop {
using namespace Eigen;
class QuadCasadiSystem : public CasadiSystem<> {
private:
  Rn<> state_manifold_;

public:
  QuadCasadiSystem(VectorXd parameters, bool use_code_generation = false);
  cs::Function casadi_step();
  cs::Function computeBodyZAxes();
};
}

#endif // QUADCASADISYSTEM_H
