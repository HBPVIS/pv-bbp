#ifndef __vtkSMBBPVolumeDensityReaderInstantiator_h
#define __vtkSMBBPVolumeDensityReaderInstantiator_h

#include "vtkInstantiator.h"



class VTK_EXPORT vtkSMBBPVolumeDensityReaderInstantiator
{
  public:
  vtkSMBBPVolumeDensityReaderInstantiator();
  ~vtkSMBBPVolumeDensityReaderInstantiator();
  private:
  static void ClassInitialize();
  static void ClassFinalize();
  static unsigned int Count;
}; 

static vtkSMBBPVolumeDensityReaderInstantiator vtkSMBBPVolumeDensityReaderInstantiatorInitializer;

#endif
