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

KeyFrame0001 = CameraKeyFrame( FocalPathPoints=[250.0, 500.0, 250.0], FocalPoint=[-2203.023551710879, 500.0, 250.0], PositionPathPoints=[2529.03, 500.0, -1000.0, 2093.7741757975127, 1839.5799825177303, -875.0, 954.2595673490897, 2667.486148213472, -750.0, -454.2581507021166, 2667.486608509874, -625.0, -1593.7733002613156, 1839.581187589426, -500.0, -2029.0299999995132, 500.0014895507461, -375.0, -1593.7750513329224, -839.5787774454627, -250.0, -454.2609839957628, -1667.485687916144, -125.0, 954.2567340548422, -1667.4870688053502, 0.0, 2093.7724247243304, -839.5823926605503, 125.0, 2529.029999998053, 499.9970208985078, 250.0, 2093.7759268675445, 1839.5775723726224, 375.0, 954.2624006421328, 2667.485227617891, 500.0, -454.25531740726706, 2667.4875290999003, 625.0, -1593.7715491865574, 1839.583597731102, 750.0, -2029.0299999956192, 500.00446865224035, 875.0, -1593.7768024013792, -839.5763672992096, 1000.0, -454.2638172882057, -1667.484767318711, 1125.0, 954.2539007593891, -1667.487989393525, 1250.0, 2093.7706736479986, -839.5848028010794, 1375.0], ClosedPositionPath=1, ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0] )

KeyFrame0002 = CameraKeyFrame( ParallelScale=1224.744871391589, Position=[2529.0272558580014, 500.0, 250.0], ViewUp=[0.0, 0.0, 1.0], KeyTime=1.0, FocalPoint=[-2203.023551710879, 500.0, 250.0] )
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
