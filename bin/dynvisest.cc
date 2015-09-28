#include <cstdio>
#include <iostream>
#include "dynvisinsview.h"
#include "viewer.h"
#include "params.h"

using namespace std;
using namespace gcop;
using namespace Eigen;

char **g_argv = 0;

Params params;

void solver_process(Viewer* viewer)
{
  viewer->animDelay = 1000;
  DynVisIns vi;

  params.GetBool("useImu", vi.useImu);
  params.GetBool("useCam", vi.useCam);
  params.GetBool("useDyn", vi.useDyn);
  params.GetBool("usePrior", vi.usePrior);

  params.GetBool("optBias", vi.optBias);

  params.GetBool("sphMeas", vi.sphMeas);

  params.GetDouble("pxStd", vi.pxStd);
  params.GetDouble("dwStd", vi.dwStd);
  params.GetDouble("dvStd", vi.dvStd);
  params.GetDouble("wStd", vi.wStd);
  params.GetDouble("aStd", vi.aStd);

  params.GetVector3d("g0", vi.g0);

  VectorXd c0(12);
  if (params.GetVectorXd("x0", c0))
    DynVisIns::ToState(vi.x0, c0);

  VectorXd P0d(12);
  if (params.GetVectorXd("P0", P0d))
    vi.x0.P.diagonal() = P0d;

  DynVisInsView view(vi);
  viewer->Add(view);

  bool sim = false;
  params.GetBool("sim", sim);

  int ni = 2;
  params.GetInt("ni", ni);

  params.GetInt("maxCams", vi.maxCams);


  double sim_dt = .125;
  double sim_tf = 4;
  params.GetDouble("sim_dt", sim_dt);
  params.GetDouble("sim_tf", sim_tf);

  DynVisIns tvi;

  DynVisInsView tview(tvi);
  tview.bodyView.rgba[0] = 1;
  tview.bodyView.rgba[1] = 0;
  tview.bodyView.rgba[2] = 0;


  if (sim) {
    //    tvi.xs = vector<Body3dState>(25);
    //    tvi.ls = vector<Vector3d>(36);
    
    VectorXd sim_c0(12);
    if (params.GetVectorXd("sim_x0", sim_c0))
      DynVisIns::ToState(tvi.x0, sim_c0);

    int ns = (int)(sim_tf/sim_dt);

    vi.SimData(tvi, ns, 25, ni, sim_dt);
    viewer->Add(tview);
  } else {
    if (!vi.LoadFile(g_argv[1])) {
      std::cerr << "ERROR: unable to open file " << g_argv[1] << "\n";
      return;
    }
  }

  getchar();
  getchar();

  vi.Compute();
  cout << "done" << endl;

  Vector12d c;
  DynVisIns::FromState(c, vi.cams[vi.camId0].x);

  cout << c.transpose() << endl;

  /*
  for (int i = 0; i< vi.xs.size(); ++i) {
    if (i < tvi.xs.size()) {
      Vector3d e;
      SO3::Instance().log(e, tvi.xs[i].R);
      cout << e.transpose() << " " << tvi.xs[i].p.transpose() << tvi.xs[i].w.transpose() << tvi.xs[i].v.transpose() << endl;
      SO3::Instance().log(e, vi.xs[i].R);
      cout << e.transpose() << " " << vi.xs[i].p.transpose() << vi.xs[i].w.transpose() << vi.xs[i].v.transpose() << endl;

    }
    //    if (vi.optBias)
    //      cout << Matrix<double, 15, 1>(vi.v + 15*i).transpose() << endl;
    //    else
    //      cout << Matrix<double, 9, 1>(vi.v + 9*i).transpose() << endl;
  }
  */

  //  for (int i = 0; i<vi.ls.size(); ++i) {
    //    cout << vi.ls[i].transpose() << endl;
    //    int i0 = vi.optBias ? vi.xs.size()*15 : vi.xs.size()*9;
    //    cout << Vector3d(vi.v + i0 + 3*i).transpose() << endl;
  //  }
  while(1)
    usleep(10000);  
}


int main(int argc, char** argv) {
  //  google::InitGoogleLogging(argv[0]);
  

  if (argc < 2) {
    std::cerr << "usage: visest input_file [options_file]\n";
    return 1;
  }

  if (argc > 2)
    params.Load(argv[2]);
  else
    params.Load("../../bin/dynvisest.cfg");

  g_argv = argv;

  //  cout << "Using params: " << vi.params.transpose() << endl;

  /*
  ofstream fout("out.txt");

  fout << vi.qs.size() << " " << vi.ps.size() << endl;
  for (int i=0; i < vi.qs.size(); ++i) {
    fout << vi.qs[i].transpose() << " ";
  }
  fout << endl;
  for (int i=0; i < vi.ps.size(); ++i) {
    fout << vi.ps[i].transpose() << " ";
  }
  fout.close();
  */

  Viewer *viewer = new Viewer;
  viewer->Init(&argc, argv);
  viewer->frameName = "body3d/frames/frame";
  viewer->SetCamera(106.5, 59, 1, -2.15, -9.3);
 
  pthread_t dummy;
  pthread_create( &dummy, NULL, (void *(*) (void *)) solver_process, viewer);

 
  viewer->Start();
  return 0;
}
