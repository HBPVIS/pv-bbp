# --script=C:\Code\plugins\pv-bbp\scripts\5K-voltage-differentials.py

from paraview.simple import *
from time import time
import os.path

#
# BBP
#
import voltageDifferentials

##########################################################
#
##########################################################
paraview.simple._DisableFirstRenderCameraReset()

filename = "snapshot"
ourvideo=[1150,1250]
config,target,report,neurons,filepath,pluginpath,frames,starttime,piston = voltageDifferentials.getOptions(filename,2,ourvideo)
RenderView1,timesteps = voltageDifferentials.initPipeline(3940,2160,config,target,report,neurons,piston)
AnimationScene1 = voltageDifferentials.helixCameraHorizontal(timesteps,frames)
voltageDifferentials.alphaKeyframe(timesteps,frames)
#voltageDifferentials.columnBox()

print "Animation time is ", AnimationScene1.AnimationTime
print "Animation end is  ", AnimationScene1.EndTime

#
# Render frames 
#
for num in range(starttime, frames[1]):
  start =  time()
  AnimationScene1.AnimationTime = timesteps.GetElement(num)
  RenderView1.ViewTime = timesteps.GetElement(num)
  Render()
  finish = time()
  timestring = "Rendered %04d %07.3f\n" % (num,finish-start)
  print timestring
  str = filepath + "/" + filename + (".%04d.png" % num)
  if (not os.path.isfile(str)):
    WriteImage(str)
  else :
    print "Skipped image ", num, ",", str

  AnimationScene1.AnimationTime = num
  RenderView1.ViewTime = num
  Render()
    
