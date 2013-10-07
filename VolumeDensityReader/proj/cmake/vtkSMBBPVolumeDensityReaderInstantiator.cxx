#include "vtkSMBBPVolumeDensityReaderInstantiator.h"
  
extern vtkObject* vtkInstantiatorvtkVolumeDensityReaderNew();


  
void vtkSMBBPVolumeDensityReaderInstantiator::ClassInitialize()
{
  
  vtkInstantiator::RegisterInstantiator("vtkVolumeDensityReader", vtkInstantiatorvtkVolumeDensityReaderNew);

  
}
          
void vtkSMBBPVolumeDensityReaderInstantiator::ClassFinalize()
{ 

  vtkInstantiator::UnRegisterInstantiator("vtkVolumeDensityReader", vtkInstantiatorvtkVolumeDensityReaderNew);

  
}

vtkSMBBPVolumeDensityReaderInstantiator::vtkSMBBPVolumeDensityReaderInstantiator()
{
  if(++vtkSMBBPVolumeDensityReaderInstantiator::Count == 1)
    { 
    vtkSMBBPVolumeDensityReaderInstantiator::ClassInitialize(); 
    }
}

vtkSMBBPVolumeDensityReaderInstantiator::~vtkSMBBPVolumeDensityReaderInstantiator()
{
  if(--vtkSMBBPVolumeDensityReaderInstantiator::Count == 0)
    { 
    vtkSMBBPVolumeDensityReaderInstantiator::ClassFinalize(); 
    }
}

// Number of translation units that include this class's header.
// Purposely not initialized.  Default is static initialization to 0.
unsigned int vtkSMBBPVolumeDensityReaderInstantiator::Count;
