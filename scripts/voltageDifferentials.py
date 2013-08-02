#!/usr/bin/python
# --script=C:\Code\plugins\pv-bbp\scripts\5K-voltage-differentials.py

from paraview.simple import *
from optparse import OptionParser
from glob import glob
import re
import socket

##########################################################
#
##########################################################
def get_pic_index(filename):
  match = re.search('\.([0-9]+).png', filename)
  try:
    index = match.group(1)
  except AttributeError:
    raise ValueError("Filename has the wrong pattern. Expected:'filename.xxxx.png', got: '{}'".format(filename))
  return int(index)
  
##########################################################
#
##########################################################
def getOptions(filename,neuronsDefault,videoMinMax):
  parser = OptionParser()
  parser.add_option("-c", "--config",  dest="config",     type="string", default="C:\\data\\bbp\\egpgv\\centralV.cfg")
  parser.add_option("-n", "--neurons", dest="neurons",    type="int",    default=neuronsDefault)
  parser.add_option("-t", "--target",  dest="target",     type="string", default="1K")
  parser.add_option("-d", "--dir",     dest="dir",        type="string", default="c:\\temp\\")
  parser.add_option("-p", "--plugin",  dest="pluginpath", type="string", default="c:\\cmakebuild\\plugin-vs10s\\bin\\Debug\\")
  parser.add_option("-s", "--start",   dest="startframe", type="int",    default=videoMinMax[0])
  parser.add_option("-e", "--end",     dest="endframe",   type="int",    default=videoMinMax[1])
  parser.add_option("-i", "--piston",  dest="piston",     type="int",    default=0)

  try :
    (options, args) = parser.parse_args()
    config     = options.config
    target     = options.target
    neurons    = options.neurons
    filepath   = options.dir
    pluginpath = options.pluginpath
    frames     = [options.startframe, options.endframe]
    piston     = options.piston
  except :
    config     = "C:\\data\\bbp\\egpgv\\centralV.cfg"
    target     = "1K"
    neurons    = 5
    filepath   = "c:\\temp\\"
    pluginpath = "c:\\cmakebuild\\plugin-vs10s\\bin\\Debug\\"
    frames     = [0,1]
    piston     = 0
    pass
    
  if (socket.gethostname()=="crusca"):
    print "Setting CRUSCA preferences"
    filepath = "D:\\temp\\"
    pluginpath = ""
  elif (socket.gethostname()=="dino"):
    print "Setting DINO preferences"
    filepath = "C:\\temp\\"
    pluginpath = ""
  else:
# the script calling this should pass options in on the command line
    servermanager.LoadPlugin(pluginpath + 'libpv_zoltan.so')
    servermanager.LoadPlugin(pluginpath + 'libpv_BBP.so')
  
  print "==============================="

  existing_indices = sorted([get_pic_index(f) for f in glob(filepath+'/'+filename+'.*.png')])
  try:
    lastindex = existing_indices[-1]
  except: 
    lastindex = frames[0]
    print "Did not find an existing file to resume from"
    
  blend_overlap = 15
  starttime = (lastindex-blend_overlap)
  if (starttime<frames[0]):
    starttime = frames[0]

  report = "voltage"  + target
  print "BlueConfig  = ", config
  print "Target      = ", target
  print "report      = ", report
  print "neurons     = ", neurons
  print "File path   = ", filepath
  print "Plugin path = ", pluginpath
  print "Start Frame = ", frames[0]
  print "End Frame   = ", frames[1]
  print "Previous N  = ", lastindex
  print "Start N     = ", starttime
  print "Piston      = ", piston
  print "==============================="
  
  return config,target,report,neurons,filepath,pluginpath,frames,starttime,piston

##########################################################
#
##########################################################
def initPipeline(xsize,ysize,config,target,report,neurons,enablepiston):
  #
  # Create BBP reader and set params
  #
  BlueConfigcircuitreader1 = BlueConfigcircuitreader()
  BlueConfigcircuitreader1.BlueConfigFileName = config
  BlueConfigcircuitreader1.DefaultTarget = target
  BlueConfigcircuitreader1.ReportName = report;
  BlueConfigcircuitreader1.MaximumNumberOfNeurons = neurons
  BlueConfigcircuitreader1.DeleteExperiment = 0
  # update information before setting point arrays and getting time steps
  BlueConfigcircuitreader1.UpdatePipelineInformation()
  BlueConfigcircuitreader1.PointArrays = ['Normal', 'RTNeuron Opacity', 'Voltage']
  timesteps = BlueConfigcircuitreader1.GetProperty("TimestepValues")

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
  DataRepresentation1.EnablePiston = enablepiston
  DataRepresentation1.Opacity = 0.9999
  DataRepresentation1.OpacityArray = 'NeuronAlpha'
  DataRepresentation1.ColorArrayName = ('POINT_DATA', 'Voltage')
  DataRepresentation1.LookupTable = a1_Voltage_PVLookupTable

  #
  # Display neurons using custom renderer and params
  #
  RenderView1 = GetRenderView()
  #RenderView1.ViewSize = [3840, 2040];
  RenderView1.ViewSize = [xsize,ysize];
  RenderView1.CenterAxesVisibility = 0
  RenderView1.Background2 = [0.0, 0.0, 0.0]
  RenderView1.CameraViewUp = [0.0, 1.0, 0.0]

  #
  # Add scalar bar to display
  #
  RenderView1.Representations.append(ScalarBarWidgetRepresentation1)

  return RenderView1,timesteps

##########################################################
#
##########################################################
def columnBox():
  #
  # Box (roughly column size)
  #
  Box1 = Box()
  Box1.Center = [250.0, 500.0, 250.0]
  Box1.ZLength = 1000.0
  Box1.YLength = 2000.0
  Box1.XLength = 1000.0

  Outline1 = OutlineCorners()
  DataRepresentation11 = Show()
  return
  
##########################################################
#
##########################################################
def helixCamera1(timesteps,trackMinMax):
  #
  # Helix camera path
  #
  AnimationScene1 = GetAnimationScene()
  AnimationScene1.Caching = False
  AnimationScene1.NumberOfFrames = len(timesteps)
  lastTime = timesteps.GetElement(len(timesteps)-1)
  AnimationScene1.EndTime = lastTime

  CameraAnimationCue1 = GetCameraTrack()
  CameraAnimationCue1.Mode = 'Path-based'

  t0 = (1.0*trackMinMax[0])/len(timesteps)
  t1 = (1.0*trackMinMax[1])/len(timesteps)
  KeyFrame0001 = CameraKeyFrame( FocalPathPoints = [262.5673476682575, 129.1111549275621, -102.35349871777954, 121.99271515964845, 750.6923459584652, 195.13225941979647, -104.45594093314287, 694.8083704276569, 155.31315835202545, 327.50085699591756, 692.3904484461164, 662.0557078675853, 507.59845183506593, 136.42294661331215, 405.5545771275172], PositionPathPoints = [2775.0650101176698, -293.8676431444712, 902.4534498544155, 1601.2444668060134, -240.42757621908976, -1232.7187383174814, -830.5054902822168, -186.9868615277025, -1385.712261091821, -2262.6793772103165, -133.54545265531937, 585.5037989115588, -1365.724213308147, -80.10394762714174, 2850.95876895307, 1027.6752118338827, -26.663033734223347, 3307.526444383927, 2695.615979389877, 26.777150510469028, 1531.354007814533, 2089.6718547201917, 80.21714396002278, -828.6558791720969, -227.63132335949902, 133.65766259143922, -1581.597360049865, -2105.034286925602, 187.09893483120837, -28.480196211010934, -1799.657397610877, 240.54048949647822, 2388.865401237491, 405.00391353263717, 293.98159321579595, 3426.3061626193103, 2462.261011589092, 347.42193125160634, 2120.738031665125, 2462.2673816268343, 400.8618996856174, -315.8199815021004, 405.0171100464495, 454.30223771791134, -1621.3988692563653, -1799.6496255525562, 507.7433414334624, -583.9696354417773, -2105.039154482476, 561.1848960982102, 1833.3743652439384, -227.64431175432978, 614.6261683411868, 3386.5013454728833, 2089.66280321252, 668.0666869765633, 2633.5719811310446, 2695.619267700825, 721.5066804271527, 273.5652624878354, 1027.6877872729522, 774.9468646689945, -1502.6158952497676, -1365.7140251002372, 828.3877785578238, -1046.0607342366402, -2262.6810344161995, 881.8292835844685, 1219.3895458479901, -830.5174544413418, 935.2706924593004, 3190.6130942582977, 1601.2333025735973, 988.711407154842, 3037.632286421661], ClosedPositionPath=1, ClosedFocalPath = 1, ParallelScale=1224.74, Position=[-3785.9285828117404, 495.54355537478784, 25.41765256980716], KeyTime=t0  )

  KeyFrame0002 = CameraKeyFrame( ParallelScale=1224.74, Position=[-3785.9285828117404, 495.54355537478784, 25.41765256980716], KeyTime=t1, FocalPoint=[250.0, 500.0, 250.0] )

  CameraAnimationCue1.KeyFrames = [ KeyFrame0001, KeyFrame0002 ]
  
  return AnimationScene1

##########################################################
#
##########################################################
def helixCameraHorizontal(timesteps,trackMinMax):
  #
  # Helix camera path
  #
  AnimationScene1 = GetAnimationScene()
  AnimationScene1.NumberOfFrames = len(timesteps)
  lastTime = timesteps.GetElement(len(timesteps)-1)
  AnimationScene1.EndTime = lastTime

  CameraAnimationCue1 = GetCameraTrack()
  CameraAnimationCue1.Mode = 'Path-based'

  t0 = (1.0*trackMinMax[0])/len(timesteps)
  t1 = (1.0*trackMinMax[1])/len(timesteps)
  KeyFrame0001 = CameraKeyFrame( FocalPathPoints=[250.0, 500.0, 250.0], FocalPoint=[-2203.023551710879, 500.0, 250.0], PositionPathPoints=[2529.03, 500.0, -1000.0, 2093.7741757975127, 1839.5799825177303, -875.0, 954.2595673490897, 2667.486148213472, -750.0, -454.2581507021166, 2667.486608509874, -625.0, -1593.7733002613156, 1839.581187589426, -500.0, -2029.0299999995132, 500.0014895507461, -375.0, -1593.7750513329224, -839.5787774454627, -250.0, -454.2609839957628, -1667.485687916144, -125.0, 954.2567340548422, -1667.4870688053502, 0.0, 2093.7724247243304, -839.5823926605503, 125.0, 2529.029999998053, 499.9970208985078, 250.0, 2093.7759268675445, 1839.5775723726224, 375.0, 954.2624006421328, 2667.485227617891, 500.0, -454.25531740726706, 2667.4875290999003, 625.0, -1593.7715491865574, 1839.583597731102, 750.0, -2029.0299999956192, 500.00446865224035, 875.0, -1593.7768024013792, -839.5763672992096, 1000.0, -454.2638172882057, -1667.484767318711, 1125.0, 954.2539007593891, -1667.487989393525, 1250.0, 2093.7706736479986, -839.5848028010794, 1375.0], ClosedPositionPath=1, ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0], KeyTime=t0 )

  KeyFrame0002 = CameraKeyFrame( ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0], KeyTime=t1, FocalPoint=[-2203.023551710879, 500.0, 250.0] )
  CameraAnimationCue1.KeyFrames = [ KeyFrame0001, KeyFrame0002 ]
  return AnimationScene1

  
##########################################################
#
##########################################################
def alphaKeyframe(timesteps,trackMinMax):
  #
  # Add a blend animation track
  #
  # if not rendering full time sequence,
  # convert time frame to local segment 
  t0 = (1.0*trackMinMax[0])/len(timesteps)
  t2 = (1.0*trackMinMax[1])/len(timesteps)
  t1 = t0+(t2-t0)/2.0
  print "keytimes are ", t0, ", ", t1, ", ", t2
  KeyFrameAnimationCue1 = GetAnimationTrack( 'DifferentialBlendFactor' )
  KeyFrame0010 = CompositeKeyFrame( KeyTime=t0,  KeyValues=[0.65] )
  KeyFrame0011 = CompositeKeyFrame( KeyTime=t1,  KeyValues=[0.65] )
  KeyFrame0012 = CompositeKeyFrame( KeyTime=t2,  KeyValues=[0.00] )
  KeyFrameAnimationCue1.KeyFrames = [ KeyFrame0010, KeyFrame0011, KeyFrame0012 ]
  return
