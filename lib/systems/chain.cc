#include "chain.h"

using namespace gcop;
using namespace Eigen;

Chain::Chain(int nb, bool fixed) : 
	Mbs(nb, nb - 1 + 6*(!fixed), fixed) {

		if (fixed)
			basetype = FIXEDBASE;

		// no gravity
		//ag << 0, 0, -9.81; 

		// structural properties
		Vector3d ds(.3, .1, .1);
		double m = 1;  
		double l0 = .3;
		double l1 = .3;
		double l2 = .3;
		Vector6d dummyvec;

		//Specifying the bodies and Joints 

		//Body0 fixed body:
		this->links[0].ds = Vector3d(0,0,0);
		this->links[0].m = 0.0001;
		dummyvec<<1e-8,1e-8,1e-8,0.0001,0.0001,0.0001;
		this->links[0].I = dummyvec;
		this->pis[0] = -1;

		//Joint 0to1
		this->joints[0].gp.setIdentity();
		//this->se3.rpyxyz2g(this->joints[0].gp, Vector3d(0,0,0), Vector3d(0.0, 0, 0.0));
		this->se3.rpyxyz2g(this->joints[0].gc, Vector3d(0,0,0), Vector3d(-0.05, 0, 0));
		this->joints[0].a.setZero();
		this->joints[0].a[2] = 1;
		this->joints[0].lower = -11.7;
		this->joints[0].upper = 11.7;    

		X.lb.r[0] = this->joints[0].lower;
		X.ub.r[0] = this->joints[0].upper;      
		X.lb.dr[0] = -15;
		X.ub.dr[0] = 15;

		U.lb[0] = -100;
		U.ub[0] = 100;


		//Body1 left_upper_shoulder:
		this->links[1].ds = Vector3d(0.1,0.05,0.05);
		this->links[1].m = 5.7044;
		dummyvec<<0.047,0.037,0.036,5.7044,5.7044,5.7044;
		this->links[1].I = dummyvec;
		this->pis[1] = 0;

		//Joint 1to2
		this->se3.rpyxyz2g(this->joints[1].gp, Vector3d(0,0,0), Vector3d(0.05, 0.0, 0.0));
		this->se3.rpyxyz2g(this->joints[1].gc, Vector3d(0,0,0), Vector3d(-0.05, 0.0, 0));
		this->joints[1].a.setZero();
		this->joints[1].a[1] = 1;
		this->joints[1].lower = -10.147;
		this->joints[1].upper = 11.047;    

		X.lb.r[1] = this->joints[1].lower;
		X.ub.r[1] = this->joints[1].upper;      
		X.lb.dr[1] = -15;
		X.ub.dr[1] = 15;

		U.lb[1] = -100;
		U.ub[1] = 100;


		//Body2 left_upper_elbow:
		this->links[2].ds = Vector3d(0.1,0.05,0.05);
		this->links[2].m = 4.313;
		dummyvec<<0.027,0.028,0.012,4.313,4.313,4.313;
		this->links[2].I = dummyvec;
		this->pis[2] = 1;

		//Joint 2to3
		this->se3.rpyxyz2g(this->joints[2].gp, Vector3d(0,0,0), Vector3d(0.05, 0.0, 0.0));
		this->se3.rpyxyz2g(this->joints[2].gc, Vector3d(0,0,0), Vector3d(-0.05, 0.0, 0));
		this->joints[2].a.setZero();
		this->joints[2].a[1] = 1;
		this->joints[2].lower = -12.618;
		this->joints[2].upper = 12.618;    

		X.lb.r[2] = this->joints[2].lower;
		X.ub.r[2] = this->joints[2].upper;      
		X.lb.dr[2] = -1.5;
		X.ub.dr[2] = 1.5;

		U.lb[2] = -100;
		U.ub[2] = 100;


		//Body3 left_lower_elbow:
		this->links[3].ds = Vector3d(0.1,0.05,0.05);
		this->links[3].m = 2.073;
		dummyvec<<0.007,0.0131,0.009,2.073,2.073,2.073;
		this->links[3].I = dummyvec;
		this->pis[3] = 2;

		U.bnd = true;
		//this->Init();
	}
