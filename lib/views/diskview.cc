#include "diskview.h"

using namespace gcop;

DiskView::DiskView(const Disk &disk, Matrix4d *g) : 
  Geom3dView("Disk", g), disk(disk)
{
  qobj = gluNewQuadric();
}


DiskView::~DiskView()
{
  free(qobj);
}

void DiskView::RenderGeom()
{  
  glPushMatrix();
  glTranslated(disk.o[0], disk.o[1], 0);
  gluDisk(qobj, 0, disk.r, 10, 10);
  glPopMatrix();
}
