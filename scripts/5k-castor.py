# --script=C:\Code\plugins\pv-bbp\scripts\5K-voltage-differentials.py

from paraview.simple import *

from optparse import OptionParser
from time import time
from glob import glob
import os.path
import re
import socket
import math
import sys

#############################
def get_pic_index(filename):
    match = re.search('\.([0-9]+).png', filename)
    try:
        index = match.group(1)
    except AttributeError:
        raise ValueError("Filename has the wrong pattern. Expected:'filename.xxxx.png', got: '{}'".format(filename))
    return int(index)
#############################

paraview.simple._DisableFirstRenderCameraReset()

parser = OptionParser()
parser.add_option("-c", "--config",  dest="config",     type="string", default="C:\\data\\bbp\\egpgv\\centralV.cfg")
parser.add_option("-n", "--neurons", dest="neurons", type="int",    default=1)
parser.add_option("-t", "--target",  dest="target",  type="string", default="1K")
parser.add_option("-d", "--dir",     dest="dir",        type="string", default="/users/biddisco/todi/")
parser.add_option("-p", "--plugin",  dest="pluginpath", type="string", default="/project/csvis/biddisco/todi/build/plugins/bin/")
parser.add_option("-s", "--start",   dest="startframe", type="int",    default=50)
parser.add_option("-e", "--end",     dest="endframe",   type="int",    default=250)

try :
  (options, args) = parser.parse_args()
  config     = options.config
  target     = options.target
  neurons    = options.neurons
  filepath   = options.dir
  pluginpath = options.pluginpath
  startframe = options.startframe
  endframe   = options.endframe
except :
  config     = "C:\\data\\bbp\\egpgv\\centralV.cfg"
  target = "1K"
  neurons = 5
  startframe = 0
  endframe   = 1
  pass

print "==============================="
if (socket.gethostname()=="crusca"):
  print "Setting CRUSCA preferences"
  filepath = "D:\\temp\\"
  pluginpath = ""
elif (socket.gethostname()=="dino"):
  print "Setting DINO preferences"
  filepath = "C:\\temp\\"
  pluginpath = ""
else:
# the script calling this passes options in on the command line
  servermanager.LoadPlugin(pluginpath + 'libpv_zoltan.so')
  servermanager.LoadPlugin(pluginpath + 'libpv_BBP.so')

filename = "snapshot"
existing_indices = sorted([get_pic_index(f) for f in glob(filepath+'/'+filename+'.*.png')])
try:
  lastindex = existing_indices[-1]
except: 
  lastindex = startframe
  print "Did not find an existing file to resume from"
  
report = "voltage"  + target
print "BlueConfig  = ", config
print "Target    = ", target
print "neurons   = ", neurons
print "report    = ", report
print "File path = ", filepath
print "Plugin path = ", pluginpath
print "Start Frame = ", startframe
print "End Frame   = ", endframe
print "Previous N  = ", lastindex
print "==============================="
  
paraview.simple._DisableFirstRenderCameraReset()
 
#
# Create BBP reader and set params
#
BlueConfigcircuitreader1 = BlueConfigcircuitreader()
BlueConfigcircuitreader1.BlueConfigFileName = config
BlueConfigcircuitreader1.DefaultTarget = target
BlueConfigcircuitreader1.ReportName = report;
BlueConfigcircuitreader1.MaximumNumberOfNeurons = neurons
BlueConfigcircuitreader1.DeleteExperiment = 0
# update information before setting point arrays 
BlueConfigcircuitreader1.UpdatePipelineInformation()
BlueConfigcircuitreader1.PointArrays = ['Normal', 'RTNeuron Opacity', 'Voltage']

#
# Colour table for scalar RTNeuron_Opacity
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
# Show the current time
#
AnnotateTimeFilter1 = AnnotateTimeFilter(BlueConfigcircuitreader1)
AnnotateTimeFilter1.Format = 'Time: %06.1f'
DataRepresentation2 = Show()
DataRepresentation2.FontSize = 12
DataRepresentation2.FontFamily = 'Courier'
DataRepresentation2.WindowLocation = 'LowerRightCorner'
DataRepresentation2.Justification = 'Center'

#
# Temporal Difference filter
#
TemporalDifferenceFilter1 = TemporalDifferenceFilter(BlueConfigcircuitreader1)
TemporalDifferenceFilter1.Arraystoprocess = []
TemporalDifferenceFilter1.Arraystoprocess = ['Voltage']
TemporalDifferenceFilter1.ComputeMagnitudes = 1

#
# Exponential Decay filter
#
ExponentialDecayFilter1 = ExponentialDecayFilter(TemporalDifferenceFilter1)
ExponentialDecayFilter1.Arraystoprocess = ['delta_Voltage']

#
# Neuron Alpha filter
#
NeuronAlphaFunction1 = NeuronAlphaFunction(ExponentialDecayFilter1)
NeuronAlphaFunction1.DifferentialVoltageArray = 'exp_delta_Voltage'
NeuronAlphaFunction1.DifferentialBlendFactor = 0.65

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
RenderView1.ViewSize = [3840, 2040];
#RenderView1.ViewSize = [1024, 768];
RenderView1.CenterAxesVisibility = 0
RenderView1.Background2 = [0.0, 0.0, 0.0]
RenderView1.CameraViewUp = [0.0, 1.0, 0.0]

#
# Add scalar bar to display
#
RenderView1.Representations.append(ScalarBarWidgetRepresentation1)

#
# Render Scene
#
Render()

#
# Helix camera path
#
timesteps = BlueConfigcircuitreader1.GetProperty("TimestepValues")
AnimationScene1 = GetAnimationScene()
AnimationScene1.NumberOfFrames = len(timesteps)
lastTime = timesteps.GetElement(len(timesteps)-1)
AnimationScene1.EndTime = lastTime

CameraAnimationCue1 = GetCameraTrack()
CameraAnimationCue1.Mode = 'Path-based'

KeyFrame0001 = CameraKeyFrame( FocalPathPoints=[250.0, 500.0, 250.0], FocalPoint=[-2203.023551710879, 500.0, 250.0], PositionPathPoints=[2529.03, 500.0, -1000.0, 2093.7741757975127, 1839.5799825177303, -875.0, 954.2595673490897, 2667.486148213472, -750.0, -454.2581507021166, 2667.486608509874, -625.0, -1593.7733002613156, 1839.581187589426, -500.0, -2029.0299999995132, 500.0014895507461, -375.0, -1593.7750513329224, -839.5787774454627, -250.0, -454.2609839957628, -1667.485687916144, -125.0, 954.2567340548422, -1667.4870688053502, 0.0, 2093.7724247243304, -839.5823926605503, 125.0, 2529.029999998053, 499.9970208985078, 250.0, 2093.7759268675445, 1839.5775723726224, 375.0, 954.2624006421328, 2667.485227617891, 500.0, -454.25531740726706, 2667.4875290999003, 625.0, -1593.7715491865574, 1839.583597731102, 750.0, -2029.0299999956192, 500.00446865224035, 875.0, -1593.7768024013792, -839.5763672992096, 1000.0, -454.2638172882057, -1667.484767318711, 1125.0, 954.2539007593891, -1667.487989393525, 1250.0, 2093.7706736479986, -839.5848028010794, 1375.0], ClosedPositionPath=1, ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0] )

KeyFrame0002 = CameraKeyFrame( ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0], KeyTime=1.0, FocalPoint=[-2203.023551710879, 500.0, 250.0] )
CameraAnimationCue1.KeyFrames = [ KeyFrame0001, KeyFrame0002 ]

#
# Add a blend animation track
#
KeyFrameAnimationCue1 = GetAnimationTrack( 'DifferentialBlendFactor' )
KeyFrame0010 = CompositeKeyFrame( KeyTime=0.0, KeyValues=[0.65] )
KeyFrame0011 = CompositeKeyFrame( KeyTime=0.5, KeyValues=[0.65] )
KeyFrame0012 = CompositeKeyFrame( KeyTime=1.0, KeyValues=[0.00] )
KeyFrameAnimationCue1.KeyFrames = [ KeyFrame0010, KeyFrame0011, KeyFrame0012 ]

#
# set animation start to required position
#
starttime = (lastindex-25)
if (starttime<startframe):
  starttime = startframe

print "Animation time is ", AnimationScene1.AnimationTime
print "Animation end is  ", AnimationScene1.EndTime

#
# Render frames 
#
timingfile=open(filepath + "/" + "timing.txt", "w+")
for num in range(starttime, endframe):
  start =  time()
  AnimationScene1.AnimationTime = timesteps.GetElement(num)
  RenderView1.ViewTime = timesteps.GetElement(num)
  Render()
  finish = time()
  timestring = "Render %08f\n" % (finish-start)
  print timestring
  timingfile.write(timestring)
  str = filepath + "/" + filename + (".%04d.png" % num)
  if (not os.path.isfile(str)):
    WriteImage(str)
  else :
    print "Skipped write of image ", str
