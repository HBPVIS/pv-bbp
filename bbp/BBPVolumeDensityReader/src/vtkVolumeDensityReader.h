/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkVolumeDensityReader.cxx,v $

  Copyright (c) Sebastien Lasserre, Blue Brain Project
  All rights reserved.

=========================================================================*/

#ifndef __vtkVolumeDensityReader_h
#define __vtkVolumeDensityReader_h 
 
#include "vtkImageAlgorithm.h"

class VTK_EXPORT vtkVolumeDensityReader : public vtkImageAlgorithm 
{
public:
    static vtkVolumeDensityReader *New();
    //vtkTypeMacro(vtkVolumeDensityReader,vtkImageAlgorithm);
    vtkTypeRevisionMacro(vtkVolumeDensityReader,vtkImageAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    //void SetSegmentVoxelMapper(float,float,float,float,float,float,float);

    //void SetWholeExtent(const long unsigned int&, int, 
    //                    const long unsigned int&, int,
    //                    const long unsigned int&, int);


    // Description:
    // Get the WholeExtent 
    vtkGetVector6Macro(WholeExtent, int);

    // Description:
    // Specify file name of the experiment file.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    vtkSetStringMacro(Target);
    vtkGetStringMacro(Target);
    
    vtkSetStringMacro(SectionType);
    vtkGetStringMacro(SectionType);
 
protected:
    vtkVolumeDensityReader();
    ~vtkVolumeDensityReader();

    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual void ExecuteData(vtkDataObject *data);

   int WholeExtent[6]; 

    float cx;
    float cy;
    float cz;
    float dimx;
    float dimy;
    float dimz;
    float res;

    char *FileName;
    char *Target;
    char *SectionType;
 
private:
    vtkVolumeDensityReader(const vtkVolumeDensityReader&); 
    void operator=(const vtkVolumeDensityReader&);

};
 
#endif

