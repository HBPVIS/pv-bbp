try: paraview.simple
except: from paraview.simple import *

from optparse import OptionParser
from time import time
import socket
import math
import sys

paraview.simple._DisableFirstRenderCameraReset()

parser = OptionParser()
parser.add_option("-n", "--neurons", dest="neurons", type="int",    default=1)
parser.add_option("-t", "--target",  dest="target",  type="string", default="1K")
parser.add_option("-p", "--path",    dest="path",    type="string", default="/users/biddisco/todi/")
(options, args) = parser.parse_args()

target   = options.target
neurons  = options.neurons
filepath = options.path

if (socket.gethostname()=="crusca"):
  filepath = "D:\\temp\\"
else:
  servermanager.LoadPlugin('/project/csvis/biddisco/todi/build/plugins/bin/libpv_zoltan.so')
  servermanager.LoadPlugin('/project/csvis/biddisco/todi/build/plugins/bin/libpv_BBP.so')

paraview.simple._DisableFirstRenderCameraReset()
 
#
# Create BBP reader and set params
#
BlueConfigcircuitreader1 = BlueConfigcircuitreader()
BlueConfigcircuitreader1.PointArrays = ['Normal', 'RTNeuron Opacity', 'Voltage']
BlueConfigcircuitreader1.DeleteExperiment = 0
BlueConfigcircuitreader1.MaximumNumberOfNeurons = neurons
BlueConfigcircuitreader1.DefaultTarget = target

#
# Colour table for Neurons
#
a1_RTNeuronOpacity_PVLookupTable = GetLookupTableForArray( "RTNeuron Opacity", 1, NanColor=[0.0, 0.0, 0.0], RGBPoints=[0.004724141, 0.0, 0.0, 1.0, 0.3497371, 1.0, 0.0, 0.0], ColorSpace='HSV', LockScalarRange=1 )

#
# Colour table for Voltage
#
a1_Voltage_PVLookupTable = GetLookupTableForArray( "Voltage", 1, RGBPoints=[-85.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039215686275, 0.498039215686275, 0.498039215686275], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1, LockScalarRange=1 )

a1_Voltage_PiecewiseFunction = CreatePiecewiseFunction( Points=[-85.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0] )

# 
# Scalar bar showing Voltage
#
ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='Voltage', Position2=[0.11161764705882349, 0.575268817204301], Enabled=1, LabelFontSize=12, LookupTable=a1_Voltage_PVLookupTable, TitleFontSize=12, Position=[0.8883823529411765, 0.17473118279569894] )

a1_Voltage_PVLookupTable.ScalarOpacityFunction = a1_Voltage_PiecewiseFunction

#
# Temporal Difference filter
#
TemporalDifferenceFilter1 = TemporalDifferenceFilter()
TemporalDifferenceFilter1.Arraystoprocess = ['Voltage']
TemporalDifferenceFilter1.ComputeMagnitudes = 1

#
# Exponential Decay filter
#
ExponentialDecayFilter1 = ExponentialDecayFilter()
ExponentialDecayFilter1.Arraystoprocess = []
ExponentialDecayFilter1.Arraystoprocess = ['delta_Voltage']

#
# Neuron Alpha filter
#
NeuronAlphaFunction1 = NeuronAlphaFunction()
NeuronAlphaFunction1.DifferentialVoltageArray = 'exp_delta_Voltage'
NeuronAlphaFunction1.DifferentialBlendFactor = 0.75

#
# Setup custom renderer and params
#
DataRepresentation1 = Show()
DataRepresentation1.Representation = 'Depth Sort Polygons'
DataRepresentation1.EnablePiston = 1
DataRepresentation1.Opacity = 0.9999
DataRepresentation1.OpacityArray = 'NeuronAlpha'
DataRepresentation1.ColorArrayName = ('POINT_DATA', 'Voltage')
DataRepresentation1.LookupTable = a1_Voltage_PVLookupTable

#
# Display neurons using custom renderer and params
#
RenderView1 = GetRenderView()
#RenderView1.ViewSize = [1920, 1200];
RenderView1.CenterAxesVisibility = 0
RenderView1.Background2 = [0.0, 0.0, 0.0]
RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
RenderView1.CameraPosition = [282.784, 703.662, 1737.1]
RenderView1.CameraClippingRange = [6.1363304921875, 6136.3304921875]
RenderView1.CameraFocalPoint = [282.784, 703.662, 1.0]
RenderView1.CameraViewAngle = 70.0
RenderView1.CameraParallelScale = 1.0
RenderView1.CenterOfRotation = [-52.8352661132812, 338.236206054688, -51.3706665039062]

#
# Set time step (1750) (t=175.0)
#
RenderView1.ViewTime = 175.0

#
# Add scalar bar to display
#
RenderView1.Representations.append(ScalarBarWidgetRepresentation1)

#
# Render Scene
#
Render()

#
# Render 30 frames for test purposes
#
timingfile=open(filepath + "timing-full.txt", "w+")
for num in range(0, 30):
  start =  time()
  Render()
  finish = time()
  timestring = "Render %08f\n" % (finish-start)
  print timestring
  timingfile.write(timestring)
  str = filepath +  ("full-%04d.png" % num)
  WriteImage(str)
