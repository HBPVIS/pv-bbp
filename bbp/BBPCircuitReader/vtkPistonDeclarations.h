#include "vtkOpenGLExtensionManager.h"
#include "vtkPistonDataObject.h"
#include "vtkPistonScalarsColors.h"

struct cudaGraphicsResource; //keeps vtkpiston namespace decl from claiming it

//----------------------------------------------------------------------------
namespace vtkpiston {
  // Forward declarations of methods defined in the cuda implementation
  void CudaGLInit();
  void CopyFromGPU(vtkPistonDataObject *id, vtkPolyData *pd);
  int QueryNumVerts(vtkPistonDataObject *id);
  int QueryVertsPer(vtkPistonDataObject *id);
  void CudaRegisterBuffer(struct cudaGraphicsResource **vboResource,
                          GLuint vboBuffer);
  void CudaTransferToGL(vtkPistonDataObject *id, unsigned long dataObjectMTimeCache,
                        vtkPistonScalarsColors *psc,
                        struct cudaGraphicsResource **vboResources,
                        bool &hasNormals, bool &hasColors);
  bool AlmostEqualRelativeAndAbs(float A, float B,
                                 float maxDiff, float maxRelDiff);

  struct PistonGLRAII
    {
    PistonGLRAII(GLbitfield mask)
      {
      glPushAttrib(mask);
      }

    ~PistonGLRAII()
      {
      glPopAttrib();
      }
    };
}

