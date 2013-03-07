/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkNeuronAlphaFunction.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef __vtkNeuronAlphaFunction_h
#define __vtkNeuronAlphaFunction_h

#include "vtkPointSetAlgorithm.h"

/*
When mapping voltage to opacity you can choose 2 approaches:
- Make hyperpolarized (voltage near -85 mV) parts transparent and depolarized (voltage around -50 mV) 
  opaque with a linear interpolation between the end points
- Make parts at resting potential (around -65 mV) transparent and interpolate until opaque as you 
  approach -85 mV and -50 mV. Ideally this should be neuron specific since each neuron has a 
  different resting potential, but we haven't worked along this line.

Additionally you can modulate the final alpha channel by the width of the morphology at that point. 
That makes thick dendrites and somas outstand against axons.

Another technique I use for the axons is to show action potentials
(spikes) instead of voltages. First you have to compute the axonal delay at each vertex 
(I don't remember the propagation speed I'm using but I can check it if you want. 
Again this should be neuron specific, but the simulation uses the same value for all, 
at least as Srikanth told me), then you have to take the spike times from the .out simulation 
file and within a reasonable large time window (propagation velocity * maximum axonal distance) 
select the active spikes and for each vertex check if there's a spike traversing it. 
I assume that spikes have some sort of "tail" with a user configurable durantion (that must be added to the
window) because at faster playbacks you need longer tails to see the spikes.
When a spike is not traversing an axon I make it fully transparent so it doesn't add clutter.
*/

class VTK_EXPORT vtkNeuronAlphaFunction : public vtkPointSetAlgorithm
{
public:
  static vtkNeuronAlphaFunction *New();
  vtkTypeMacro(vtkNeuronAlphaFunction,vtkPointSetAlgorithm);

  // Description:
  // 
  vtkSetMacro(DifferentialBlendFactor,double);
  vtkGetMacro(DifferentialBlendFactor,double);

  // Description:
  // hyperpolarized (voltage near -85 mV)
  vtkSetMacro(HyperPolarizedVoltage,double);
  vtkGetMacro(HyperPolarizedVoltage,double);

  // Description:
  // depolarized (voltage near -50 mV)
  vtkSetMacro(DePolarizedVoltage,double);
  vtkGetMacro(DePolarizedVoltage,double);

  // Description:
  // resting potential (voltage near -65 mV)
  vtkSetMacro(RestingPotentialVoltage,double);
  vtkGetMacro(RestingPotentialVoltage,double);
  
  // Description:
  // Peak deviation/differential voltage reference
  // which will map to full opacity
  vtkSetMacro(PeakDifferentialVoltage,double);
  vtkGetMacro(PeakDifferentialVoltage,double);
  
  // Description:
  // depolarized (voltage near -50 mV)
  vtkSetMacro(VoltageTransparencyMode,int);
  vtkGetMacro(VoltageTransparencyMode,int);

  vtkSetStringMacro(Array1Name);
  vtkGetStringMacro(Array1Name);

  vtkSetStringMacro(Array2Name);
  vtkGetStringMacro(Array2Name);

  vtkSetStringMacro(Array3Name);
  vtkGetStringMacro(Array3Name);

  vtkSetStringMacro(Array4Name);
  vtkGetStringMacro(Array4Name);

protected:
   vtkNeuronAlphaFunction();
  ~vtkNeuronAlphaFunction();
  //
  double HyperPolarizedVoltage;
  double DePolarizedVoltage;
  double RestingPotentialVoltage;
  double PeakDifferentialVoltage;
  int    VoltageTransparencyMode;
  double DifferentialBlendFactor;
  //
  char  *Array1Name;
  char  *Array2Name;
  char  *Array3Name;
  char  *Array4Name;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkNeuronAlphaFunction(const vtkNeuronAlphaFunction&);  // Not implemented.
  void operator=(const vtkNeuronAlphaFunction&);  // Not implemented.
};

#endif
