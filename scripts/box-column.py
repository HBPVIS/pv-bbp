# --script=C:\Code\plugins\pv-bbp\scripts\box-column.py
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

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
parser.add_option("-n", "--neurons", dest="neurons", type="int",    default=1)
parser.add_option("-t", "--target",  dest="target",  type="string", default="1K")
parser.add_option("-p", "--path",    dest="path",    type="string", default="/users/biddisco/todi/")
try :
  (options, args) = parser.parse_args()
  target   = options.target
  neurons  = options.neurons
  filepath = options.path
except :
  target = "1K"
  neurons = 5
  pass

print "==============================="
if (socket.gethostname()=="crusca"):
  filepath = "D:\\temp\\"
elif (socket.gethostname()=="dino"):
  print "Setting DINO preferences"
  filepath = "C:\\temp\\"
#  servermanager.LoadPlugin('pv_zoltan.dll')
#  servermanager.LoadPlugin('pv_BBP.dll')
else:
  filepath="/scratch/castor/biddisco/temp/"
  servermanager.LoadPlugin('/project/csvis/biddisco/todi/build/plugins/bin/libpv_zoltan.so')
  servermanager.LoadPlugin('/project/csvis/biddisco/todi/build/plugins/bin/libpv_BBP.so')

filename = "snapshot"
existing_indices = sorted([get_pic_index(f) for f in glob(filepath+'/'+filename+'.*.png')])
try:
  lastindex = existing_indices[-1]
except: 
  lastindex = 0
  print "Did not find an existing file to resume from"
  
#
# Box (roughly column size)
#

Box1 = Box()
Box1.ZLength = 1000.0
Box1.YLength = 2000.0
Box1.XLength = 1000.0
Box1.Center = [250.0, 500.0, 250.0]

DataRepresentation1 = Show()
DataRepresentation1.Representation = 'Wireframe'

#
# Show the current time
#
AnnotateTimeFilter1 = AnnotateTimeFilter()
AnnotateTimeFilter1.Format = 'Time: %06.1f'
DataRepresentation2 = Show()
DataRepresentation2.FontSize = 12
DataRepresentation2.FontFamily = 'Courier'
DataRepresentation2.WindowLocation = 'LowerRightCorner'

RenderView1 = GetRenderView()
RenderView1.CenterOfRotation = [250.0, 500.0, 250.0]
RenderView1.CameraPosition = [250.0, 500.0, 4982.050807568878]
RenderView1.CameraFocalPoint = [250.0, 500.0, 250.0]
RenderView1.CameraClippingRange = [3689.730299493189, 6055.531569682411]
RenderView1.CameraParallelScale = 1224.744871391589


#
# Helix camera path
#
AnimationScene1 = GetAnimationScene()
AnimationScene1.NumberOfFrames = 1000

CameraAnimationCue1 = GetCameraTrack()
CameraAnimationCue1.Mode = 'Path-based'

KeyFrame0001 = CameraKeyFrame( FocalPathPoints = [262.5673476682575, 129.1111549275621, -102.35349871777954, 121.99271515964845, 750.6923459584652, 195.13225941979647, -104.45594093314287, 694.8083704276569, 155.31315835202545, 327.50085699591756, 692.3904484461164, 662.0557078675853, 507.59845183506593, 136.42294661331215, 405.5545771275172], PositionPathPoints = [2775.0650101176698, -293.8676431444712, 902.4534498544155, 1601.2444668060134, -240.42757621908976, -1232.7187383174814, -830.5054902822168, -186.9868615277025, -1385.712261091821, -2262.6793772103165, -133.54545265531937, 585.5037989115588, -1365.724213308147, -80.10394762714174, 2850.95876895307, 1027.6752118338827, -26.663033734223347, 3307.526444383927, 2695.615979389877, 26.777150510469028, 1531.354007814533, 2089.6718547201917, 80.21714396002278, -828.6558791720969, -227.63132335949902, 133.65766259143922, -1581.597360049865, -2105.034286925602, 187.09893483120837, -28.480196211010934, -1799.657397610877, 240.54048949647822, 2388.865401237491, 405.00391353263717, 293.98159321579595, 3426.3061626193103, 2462.261011589092, 347.42193125160634, 2120.738031665125, 2462.2673816268343, 400.8618996856174, -315.8199815021004, 405.0171100464495, 454.30223771791134, -1621.3988692563653, -1799.6496255525562, 507.7433414334624, -583.9696354417773, -2105.039154482476, 561.1848960982102, 1833.3743652439384, -227.64431175432978, 614.6261683411868, 3386.5013454728833, 2089.66280321252, 668.0666869765633, 2633.5719811310446, 2695.619267700825, 721.5066804271527, 273.5652624878354, 1027.6877872729522, 774.9468646689945, -1502.6158952497676, -1365.7140251002372, 828.3877785578238, -1046.0607342366402, -2262.6810344161995, 881.8292835844685, 1219.3895458479901, -830.5174544413418, 935.2706924593004, 3190.6130942582977, 1601.2333025735973, 988.711407154842, 3037.632286421661], ClosedPositionPath=1, ClosedFocalPath = 1, ParallelScale=1224.74, Position=[-3785.9285828117404, 495.54355537478784, 25.41765256980716] )

KeyFrame0002 = CameraKeyFrame( ParallelScale=1224.74, Position=[-3785.9285828117404, 495.54355537478784, 25.41765256980716], KeyTime=1.0, FocalPoint=[250.0, 500.0, 250.0] )

CameraAnimationCue1.KeyFrames = [ KeyFrame0001, KeyFrame0002 ]

#
# Display 
#
RenderView1 = GetRenderView()
RenderView1.ViewSize = [1024, 768];

Render()

#
# Render N frames for test purposes
#
timingfile=open(filepath + "/" + "timing.txt", "w+")
for num in range(0, 10):
  start =  time()
  AnimationScene1.AnimationTime = num
  RenderView1.ViewTime = num
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
