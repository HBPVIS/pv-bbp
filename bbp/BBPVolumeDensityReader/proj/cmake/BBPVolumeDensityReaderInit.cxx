#include "vtkClientServerInterpreter.h"

#ifndef PARAVIEW_BUILD_SHARED_LIBS
#define PARAVIEW_BUILD_SHARED_LIBS
#endif
#if defined(PARAVIEW_BUILD_SHARED_LIBS) && defined(_WIN32)
# define VTK_WRAP_CS_EXPORT __declspec(dllexport)
#else
# define VTK_WRAP_CS_EXPORT
#endif

#if 0
static int BBPVolumeDensityReader_NewInstance(vtkClientServerInterpreter *arlu,
        const char *type, vtkClientServerID id);

int vtkVolumeDensityReaderCommand(vtkClientServerInterpreter *, vtkObjectBase *, const char *, const vtkClientServerStream&, vtkClientServerStream& resultStrem);
vtkObjectBase *vtkVolumeDensityReaderClientServerNewCommand();

#endif

extern void vtkVolumeDensityReader_Init(vtkClientServerInterpreter* csi);


extern "C" void VTK_WRAP_CS_EXPORT BBPVolumeDensityReader_Initialize(
  vtkClientServerInterpreter *csi)
{
    vtkVolumeDensityReader_Init(csi);

#if 0
  arlu->AddNewInstanceFunction( BBPVolumeDensityReader_NewInstance);

    arlu->AddCommandFunction("vtkVolumeDensityReader",vtkVolumeDensityReaderCommand);

#endif
}

#if 0
static int BBPVolumeDensityReader_NewInstance(vtkClientServerInterpreter *arlu,
                                   const char *type, vtkClientServerID id)
{

    if (!strcmp("vtkVolumeDensityReader",type))
      {
      vtkObjectBase *ptr = vtkVolumeDensityReaderClientServerNewCommand();
      arlu->NewInstance(ptr,id);
      return 1;
      }


  return 0;
}
#endif
