#include "view.h"
#include <iostream>
#include <cstring>

using namespace gcop;

static int s_view_id = 0;

View::View(const char* name)
{
  if (name)
    strcpy(this->name, name);
  else
    this->name[0] = 0;

  this->id = (s_view_id++);

  pthread_mutex_init(&mut, 0);
}

View::~View()
{
  pthread_mutex_unlock(&mut);
  pthread_mutex_destroy(&mut);
}

void View::Lock()
{
  pthread_mutex_lock(&mut);
}

void View::Unlock()
{
  pthread_mutex_unlock(&mut);
}


void View::Render()
{
  int i=0;
  while(RenderFrame(i))
    ++i;
}


bool View::RenderFrame(int i)
{
  std::cerr << "Warning: View::RenderFrame: #" << i << " no rendering for view id#" << this->id << " :(" << name << ")" << std::endl;
  return false;
}

