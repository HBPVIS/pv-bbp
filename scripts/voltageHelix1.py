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
ourvideo=[0,5000]
config,target,report,neurons,filepath,pluginpath,startframe,endframe,starttime,piston = voltageDifferentials.getOptions(filename,2,ourvideo)
RenderView1,timesteps = voltageDifferentials.initPipeline(1024,768,config,target,report,neurons,piston)
AnimationScene1 = voltageDifferentials.helixCamera1(timesteps,ourvideo)
#AnimationScene1 = voltageDifferentials.helixCameraHorizontal(timesteps,ourvideo)
voltageDifferentials.alphaKeyframe(timesteps,ourvideo)

print "Animation time is ", AnimationScene1.AnimationTime
print "Animation end is  ", AnimationScene1.EndTime

#
# Render frames 
#
for num in range(starttime, endframe):
  AnimationScene1.AnimationTime = timesteps.GetElement(num)
  RenderView1.ViewTime = timesteps.GetElement(num)
  print "Rendering ", num
  Render()
  str = filepath + "/" + filename + (".%04d.png" % num)
  if (not os.path.isfile(str)):
    WriteImage(str)
  else :
    print "Skipped image ", num, ",", str
