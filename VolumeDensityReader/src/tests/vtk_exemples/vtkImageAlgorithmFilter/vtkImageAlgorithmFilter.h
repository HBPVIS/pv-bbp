#ifndef __vtkImageAlgorithmFilter_h
#define __vtkImageAlgorithmFilter_h
 
#include "vtkImageAlgorithm.h"
 
class VTK_IMAGING_EXPORT vtkImageAlgorithmFilter : public vtkImageAlgorithm 
{
public:
    static vtkImageAlgorithmFilter *New();
    vtkTypeRevisionMacro(vtkImageAlgorithmFilter,vtkImageAlgorithm);
    //void PrintSelf(ostream& os, vtkIndent indent);

    vtkGetVector6Macro(WholeExtent, int);
  
 
 
protected:
    vtkImageAlgorithmFilter();
    ~vtkImageAlgorithmFilter();

    int WholeExtent[6];
 
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual void ExecuteData(vtkDataObject *data);

private:
    vtkImageAlgorithmFilter(const vtkImageAlgorithmFilter&);  // Not implemented.
    void operator=(const vtkImageAlgorithmFilter&);  // Not implemented.
 
};
 
#endif
