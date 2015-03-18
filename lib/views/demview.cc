#include "demview.h"
#include <stdlib.h>
#include "GL/glut.h"
#include "utils.h"
#include "viewer.h"
#include <iostream>

using namespace gcop;
using namespace std;

DemView::DemView(const Dem& dem) :
  View("Dem"),
  dem(dem)
{
  //  color[0] = .1;
  //  color[1] = .8;
  //  color[2] = .3;
  color[0] = .6;
  color[1] = .6;
  color[2] = .6;

  wire = false;
  mesh_inds = 0;
  vertices = 0;
  normals = 0;
  texture = 0;
  Init();
}

void DemView::Reset()
{
  //  delete[] mesh_inds;
  // delete[] vertices;
  Init();
}



void DemView::Init()
{
  int i, j;
  double xyz[3];
  int ind = 0;
  //  double A[3], B[3], N[3];

  if (!vertices)
    vertices = (GLfloat*)calloc(dem.ni*dem.nj*12, sizeof(GLfloat));

  if (!normals)
    normals = (GLfloat*)calloc(dem.ni*dem.nj*3, sizeof(GLfloat));

  for (i = 0; i < dem.ni-1; ++i) {
    for (j = 0; j < dem.nj-1; ++j) {
      dem.Get(xyz, i, j);
  
      /******* normals ***********/
      double p1[3];
      double p2[3];
      p1[0] = xyz[0] + dem.cs/2;
      p1[1] = xyz[1];
      p1[2] = dem.Get(p1[0], p1[1]);

      p2[0] = xyz[0];
      p2[1] = xyz[1] + dem.cs/2;
      p2[2] = dem.Get(p2[0], p2[1]);
      
      double v1[3];
      double v2[3];
      MINUS3(v1, p1,xyz);
      MINUS3(v2, p2,xyz);
      double n[3];
      CROSS(n, v1, v2);
      double nn = NORM3(n);
      DIV3(normals + 3*(i*dem.nj +j), n, nn);
      /****************************/
    
      // The texture coordinates
      vertices[ind++] = 0; 
      vertices[ind++] = 0;
      
      //	  float c = (xyz[2] - map->min)/(map->max - map->min);
      
      /*
	vertices[ind++] = .2 + .8*c;
	vertices[ind++] = .2 + .8*c;
	vertices[ind++] = .2 + .8*c;
      */      
      vertices[ind++] = color[0];
      vertices[ind++] = color[1];
      vertices[ind++] = color[2];
      vertices[ind++] = 0;

      if(0) {
	/*
	  if (i > 0 && i < dem->ind_width - 1 && j > 0 && j < map->ind_height - 1) {
	  dem_get_ind(dem, A, i+1, j); 
	  dem_get_ind(dem, B, i, j+1); 
	  cross(N, sub(A,xyz),sub(B,xyz));
	  vertices[ind++] = N[0];
	  vertices[ind++] = N[1];
	  vertices[ind++] = N[2];
	*/
      } else {
	// The Normal vector
	vertices[ind++]=0;
	vertices[ind++]=0;
	vertices[ind++]=1;
      }
      
      // The vertex position
      vertices[ind++] = xyz[0];
      vertices[ind++] = xyz[1];
      vertices[ind++] = xyz[2];
    }
  }
  
  mesh_ind_count = (dem.nj - 1)*(dem.ni - 1)*4;
  if (!mesh_inds)
    mesh_inds = (GLuint*)calloc(mesh_ind_count, sizeof(GLuint));
  
  ind = 0;
  for (i = 0; i < dem.ni - 1; ++i) {
    for (j = 0; j < dem.nj - 1; ++j) {
      mesh_inds[ind++] = ((i+1)*dem.nj + j);
      mesh_inds[ind++] = (i*dem.nj + j);
      mesh_inds[ind++] = (i*dem.nj + j);
      mesh_inds[ind++] = (i*dem.nj + j + 1);
    }
  }
}

DemView::~DemView()
{
  delete[] mesh_inds;
  delete[] vertices;
  delete[] normals;

  if (texture)
    glDeleteTextures( 1, &texture );
}

void DemView::SetTexture(const char *fname)
{
  //  glEnable( GL_TEXTURE_2D );
  texture = Viewer::LoadTexture(fname, true);
  //  std::cout << "texture=" << texture << std::endl;
  glBindTexture( GL_TEXTURE_2D, texture );
}


bool DemView::RenderFrame(int i)
{
  /*
  glInterleavedArrays(GL_T2F_C4F_N3F_V3F, 0, vertices);
  //  glDrawElements(GL_TRIANGLES, mapFillIndexCount, GL_UNSIGNED_INT, 1);
  glDrawElements(GL_LINES, mesh_ind_count, GL_UNSIGNED_INT, mesh_inds);
  // glDrawElements(GL_TRIANGLES, mesh_ind_count, GL_UNSIGNED_INT, mesh_inds);
  */

  glColor3fv(color);
  //  glColor4f(1,1,1,0.1);  
  
  if (texture) {
    glColor4f(1,1,1,1);  
    glDisable(GL_LIGHTING);
  } else {
    glEnable(GL_LIGHTING);
  }
  //  Viewer::SetColor(color[0], color[1], color[2], 0);

  // double c = .5;
  // Viewer::SetMaterial(c,c,c, c,c,c, c,c,c,5);
  // glColor3f(.5,.5,.5);

  if (wire) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    //    glDisable(GL_DEPTH_TEST);
    //    glCullFace(GL_BACK);
  }
  double p[3];
  for (int i = 0; i < dem.ni - 1; ++i) {
        //Makes OpenGL draw a triangle at every three consecutive vertices
    glBegin(GL_TRIANGLE_STRIP);
    for (int j = 0; j < dem.nj; ++j) {
      dem.Get(p, i, j);
      
      // glNormal3fv(normals + 3*(i*dem.nj+j));
      // glNormal3dv(dem.normals + 3*(i*dem.nj+j));

      glVertex3f(p[0], p[1], p[2]);
      if (texture)
        glTexCoord2d( (p[0]-dem.o[0])/dem.w, (dem.h - (p[1]-dem.o[1]) )/dem.h);      
      
      dem.Get(p, i+1, j);

      //      glNormal3fv(normals + 3*((i+1)*dem.nj+j));
      glNormal3dv(dem.normals + 3*((i+1)*dem.nj+j));
      
      glVertex3f(p[0], p[1], p[2]);

      if (texture)
        glTexCoord2d( (p[0]-dem.o[0])/dem.w, (dem.h - (p[1]-dem.o[1]) )/dem.h);     
    }
    glEnd();
  }

  if (wire)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

      // enable blending
    glEnable(GL_BLEND);

    // enable read-only depth buffer
    glDepthMask(GL_FALSE);

    // set the blend function
    // to what we use for transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // set back to normal depth buffer mode (writable)
    glDepthMask(GL_TRUE);

    // disable blending
    glDisable(GL_BLEND);


  if (texture)
    glEnable(GL_LIGHTING);
  return false;
}

void DemView::Render()
{
  RenderFrame(0);
}

