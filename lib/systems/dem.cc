#include "dem.h"
#include <cmath>
#include "utils.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace gcop;
using namespace std;

Dem::Dem(const Dem &dem) :
  w(dem.w), h(dem.h), cs(dem.cs), ds(dem.ds), ni(dem.ni), nj(dem.nj), eps(dem.eps)
{
  memcpy(o, dem.o, 3*sizeof(double));
  data = new double[nj*ni];
  memcpy(data, dem.data, nj*ni*sizeof(double));
  odata = 0;
  normals = new double[nj*ni*3];
  memcpy(normals, dem.normals, nj*ni*3*sizeof(double));
}

Dem::Dem() : w(0), h(0), cs(0), ds(0), ni(0), nj(0), data(0), odata(0), normals(0), eps(1e-10)
{
  memset(this->o, 0, 3*sizeof(double));
}

Dem::Dem(double w, double h, double cs, double ds, const double *o) :
  w(w), h(h), cs(cs), ds(ds), eps(1e-10)
{
  if (o)
    memcpy(this->o, o, 3*sizeof(double));
  else
    memset(this->o, 0, 3*sizeof(double));
  nj = (int)(w/cs+1);
  ni = (int)(h/cs+1);
  data = new double[nj*ni];
  odata = 0;
  normals = new double[nj*ni*3];
}


Dem::Dem(const char *fname, double cs, double ds, const double *o) :
  w(0), h(0), cs(cs), ds(ds), eps(1e-10)
{
  if (o)
    memcpy(this->o, o, 3*sizeof(double));
  else
    memset(this->o, 0, 3*sizeof(double));
  Load(fname);
  odata = 0;

  ComputeNormals();
}


Dem::~Dem()
{
  delete[] normals;
  delete[] data;
  delete[] odata;
}

const double* Dem::GetNormal(int i, int j) const
{
  return normals + 3*(j + nj*i);
}

double Dem::bilinterp(const double* z, int w, int h, double xi, double yi, double eps)
{
  double t,u;
  int x1, x2, y1, y2;
  
  //  cout << "xi=" << xi << " yi=" << yi << endl;

  if (xi < 0 && xi > -eps)
    xi = 0;
  if (xi > w - 1 && xi < w - 1 + eps)
    xi = w - 1 - eps;

  if (yi < 0 && yi > -eps)
    yi = 0;
  if (yi > h - 1 && yi < h - 1 + eps )
    yi = h - 1 - eps;

  x1 = (int)floor(xi);
  y1 = (int)floor(yi);
  x2 = (int)ceil(xi);
  y2 = (int)ceil(yi);


  if (x1 < 0 || x1 >= w || y1 < 0 || y1 >= h) {     
    fprintf(stderr, "Error - dem:bilinterp - invalid (x1,y1)=(%d,%d) (w=%d, h=%d)\n", x1,y1,w, h);
    return 0;
  }
  if (x2 < 0 || x2 >= w || y2 < 0 || y2 >= h) {
    fprintf(stderr, "Error - dem:bilinterp - invalid (x2,y2)=(%d,%d) (w=%d, h=%d)\n", x2,y2, w, h);
    return 0;
  }

  t = (xi - x1);
  u = (yi - y1);

  return (1-t)*(1-u)*z[y1*w + x1] + t*(1-u)*z[y1*w + x2] + t*u*z[y2*w + x2] + (1-t)*u*z[y2*w + x1];
}



void Dem::ComputeNormals()
{
  if (!normals) {
    std::cerr << "[W] dgc::Dem::ComputeNormals: normals not initialized! Initializing now... OK" << std::endl;
    normals = new double[nj*ni*3];
  }

  double p[3]; // point
  double vx[3];  // vector in x-direction
  double vy[3];  // vector in y-direction
  double n[3];   // normal

  for (int i = 0; i < ni; ++i) {
    int j = nj-1;
    int ind = 3*(j + nj*i);
    normals[ind] = 0;
    normals[ind+1] = -1;
    normals[ind+2] = 0;
  }

  for (int j = 0; j < nj; ++j) {
    int i = 0;
    int ind = 3*(j + nj*i);
    normals[ind] = -1;
    normals[ind+1] = 0;
    normals[ind+2] = 0;
  }

  for (int i = 1; i < ni; ++i) {
    for (int j = 0; j < nj-1; ++j) {
      Get(p, i, j);
      Get(vx, i, j+1);
      Get(vy, i-1, j);
      MINUS3(vx, vx, p);
      MINUS3(vy, vy, p);
      CROSS(n, vx, vy);

      double nn = NORM3(n);

      //      cout << "vx= " << vx[0] << " " << vx[1] << " " << vx[2] << endl;
      //      cout << "vy= " << vy[0] << " " << vy[1] << " " << vy[2] << endl;
      //      cout << "n= " << n[0]/nn << " " << n[1]/nn << " " << n[2]/nn << endl;

      int ind = 3*(j + nj*i);

      normals[ind] = n[0]/nn;
      normals[ind+1] = n[1]/nn;
      normals[ind+2] = n[2]/nn;
    }
  }

  //  cout << ni << " " <<  nj << endl;
}

void Dem::Dilate(double r, bool cube) 
{
  if (!odata) {
    odata = new double[ni*nj];
  }
  memcpy(odata, data, ni*nj*sizeof(double));

  // only cube supported for now
  assert(cube);

  int di = (int)ceil(r/cs);
  if (di <= 0)
    return;

  cout << "Dem::Dilate: dilating with di=" << di << endl;

  for (int i = di; i < ni - di; ++i) {
    for (int j = di; j < nj - di; ++j) {
      // double z = odata[i*nj + j]; // height at center
      for (int k = -di; k <= di; ++k) {
        for (int l = -di; l <= di; ++l) {
          int ik = i + k;
          int jl = j + l;
          data[ik*nj + jl] = max(data[ik*nj + jl], odata[ik*nj + jl]);
        }        
      }
    }
  }
  ComputeNormals();
  //  memcpy(odata, data, ni*nj*sizeof(double));
}

void Dem::Convolve(double sigma, bool cn, double thresh)
{
  if (!odata) {
    odata = new double[ni*nj];
    memcpy(odata, data, ni*nj*sizeof(double));
  }

  int a = (int)round(2*sigma/cs);

  for (int i = 0; i < ni; ++i) {
    for (int j = 0; j < nj; ++j) {

      //      int jl = MAX(j-a,0);
      //      int jr = MIN(j+a,nj-1);
      //     int il = MAX(i-a,0);
      //     int ir = MIN(i+a,ni-1);
      
      double z = 0;
      double n = 0;
      
      //      for (int ci = il; ci <= ir; ++ci) {
      //        for (int cj = jl; cj <= jr; ++cj) {
      for (int ci = i-a; ci <= i+a; ++ci) {
        for (int cj = j-a; cj <= j+a; ++cj) {
          
          double &d = odata[cj + nj*ci];
          if (d < thresh)
            continue;
          double r =  (cj-j)*(cj-j) + (ci-i)*(ci-i);
          double g = 1/(sigma*sqrt(2*M_PI))*exp(-r*cs*cs/(2*sigma*sigma));
          n += g;
          if (ci >= 0 && ci < ni && cj >= 0 && cj < nj)
            z += g*d;
        }
      }
      data[j + nj*i] = z/n;
    }
  }
  if (cn)
    ComputeNormals();
}



bool Dem::Inside(double x, double y, double z) const
{
  if (IsValid(x,y))
    return Get(x,y) > z;
  return false;
}


void Dem::Set(int i, int j, double z)
{
  //  cout << x << " " << y << endl;
  assert(i>= 0  && i <ni);
  assert(j>= 0  && j <nj);
  data[i*nj + j] = z;
}

void Dem::Set(double x, double y, double z, double s)
{
  //  cout << x << " " << y << endl;
  if (!IsValid(x, y)) {
    std::cerr << "Warning: Dem::Set: invalid (x,y)=(" << x << "," << y << ")" << std::endl;
  }
  
  if (s < eps) {
    data[(int)((x - o[0])/cs) + nj*((int)((h - y - o[1])/cs))] = z;
  } else {
    for (double xs = x - s; xs <= x + s; xs += cs)
      for (double ys = y - s; ys <= y + s; ys += cs)
        if (IsValid(xs, ys))
          data[(int)((xs - o[0])/cs) + nj*((int)((h - ys - o[1])/cs))] = z;
  }
}

void Dem::Clear()
{
  memset(data, 0, ni*nj*sizeof(double));
}

const double* Dem::GetNormal(double x, double y) const
{
  if (!IsValid(x, y)) {
    return 0;
  }
  int ind = 3*((int)((x - o[0])/cs) + nj*((int)((h - y - o[1])/cs)));
  return normals + ind;
}


double Dem::GetNormal(double n[3], double x, double y) const
{
  if (!IsValid(x, y)) {
    n[0] = 0;
    n[1] = 0;
    n[2] = 1;
    return 0;
  }
  int ind = 3*((int)((x - o[0])/cs) + nj*((int)((h - y - o[1])/cs)));
  SET3(n, normals+ind);
  
  return o[2] + bilinterp(data, nj, ni, (x - o[0])/cs, (h - y - o[1])/cs);
}

double Dem::Get(double x, double y) const
{
  //  cout << x << " " << y << endl;
  if (!IsValid(x, y)) {
    //    std::cerr << "[W] dgc::Dem::Get: invalid (x,y)=(" << x << "," << y << ")" << std::endl;
    return 0;
  }
  return o[2] + bilinterp(data, nj, ni, (x - o[0])/cs, (h - y - o[1])/cs);
}

void Dem::Point2Index(int &i, int &j, double x, double y) const
{
  j = (x - o[0])/cs;
  i = (h - y - o[1])/cs;
}

void Dem::Index2Point(double &x, double &y, int i, int j) const
{
  x = j*cs + o[0];
  y = h - o[1] - i*cs;
}


void Dem::Get(double *p, int i, int j) const
{
  if (i < 0 || i >= ni ||
      j < 0 || j >= nj) {
    std::cerr << "Warning: Dem::Get: invalid (i,j)=(" << i << "," << j << ")" << std::endl;
    return;
  }

  p[0] = j*cs + o[0];
  p[1] = h - i*cs + o[1];
  p[2] = data[i*nj + j] + o[2];
}

void Dem::Scale(double s)
{
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      data[i*nj + j] *= s;
}


void Dem::AddBoundary(double h)
{
  for (int i = 0; i < ni; ++i) {
    data[0 + nj*i] = h;
    data[nj-1 + nj*i] = h;
  }
  for (int j = 0; j < nj; ++j) {
    data[j + nj*0] = h;
    data[j + nj*(ni-1)] = h;
  } 
}


void Dem::Load(const char *fname)
{
  int i, size;
  char* chdata;
  fstream fstr;

  fstr.open(fname, std::ios::in);
  
  if (!fstr.good()) {
    std::cerr << "[E]: Dem:Load - failed to open " << string(fname) << std::endl;
    return;
  }
  
  char line[256];
  
  fstr.getline(line, 256);
  if (line[0] != 'P' || line[1] != '6') {
    cerr << "[E]: Dem:Load - file " << fname << " should be in P6 binary format" << endl;
    return;
  }
  
  // skip comments
  do {
    fstr.getline(line, 256);
  } while(line[0] == '#');
  
  sscanf(line, "%d %d 255\n", &nj, &ni);
  
  size = nj*ni;
  data = new double[size]; // (double*)calloc(size, sizeof(double));
  normals = new double[size*3];
  chdata = new char[3*size]; // (char*)malloc(size*3);
  
  fstr.read(chdata, size*3);
  //  fread(chdata, sizeof(char), size*3, file);
  for (i = 0; i < size; i++) {
    data[i] = (unsigned char)chdata[3*i]*ds/255.0;
  }
  delete[] chdata;
  //  fclose(file);

  fstr.close();

  w = (nj - 1)*cs;
  h = (ni - 1)*cs;
}
