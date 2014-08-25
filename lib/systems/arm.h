#ifndef GCOP_ARM_H
#define GCOP_ARM_H

namespace gcop {

  /**
   * 3-dof arm (rotating base, and two links). Model assumes first joint motor is at center of mass.
   * The local frame is z-up, and arm is zero position means that is straightened along the x+ axis.
   * To account for position offset simply change the parameter p.
   */
  class Arm {
  public:
    Arm();

    virtual ~Arm();

    /**
     * Forward kinematics: 
     * @param p end-effector position
     * @param a angles: yaw, joint1, joint2
     */
    bool Fk(double p[3], const double a[3]);

    /**
     * Inverse kinematics
     * @param a angles: yaw, joint1, joint2; a[0][*] is first solution, a[1][*] is second solution
     * @param p end-effector position
     */    
    bool Ik(double a[2][3], const double p[3]);

    double l1;  ///< first link length
    double l2;  ///< second link length
    double x1;  ///< relative offset to second motor in x direction
    

  };
}

#endif

