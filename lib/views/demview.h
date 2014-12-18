#ifndef DemVIEW_H
#define DemVIEW_H

#include "view.h"
#include "dem.h"

#include "GL/glu.h"

namespace gcop {

  class DemView : public View {
  public:
    
    DemView(const Dem& dem);
    virtual ~DemView();
    
    void Reset();

    void Render();
    bool RenderFrame(int i);

    void SetTexture(const char *fname);

    float color[3];

    bool wire;          ///< draw as wire mesh (false by default)

  protected:
    void Init();
    //    GLUquadricObj *qobj;
    const Dem& dem;

    float* vertices; // GL_T2F_C4F_N3F_V3F
    int mesh_ind_count;
    GLuint* mesh_inds;

    float *normals;

    GLuint texture;     ///< texture

  };
};




#endif
