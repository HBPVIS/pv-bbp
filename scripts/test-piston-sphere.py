from paraview.simple import *
from time import time

paraview.simple._DisableFirstRenderCameraReset()

Sphere1 = Sphere()
Sphere1.ThetaResolution = 64
Sphere1.PhiResolution = 32

Elevation1 = Elevation()
Elevation1.LowPoint = [-0.499358266592026, -0.499358266592026, -0.5]
Elevation1.HighPoint = [0.499358266592026, 0.499358266592026, 0.5]

a1_Elevation_PVLookupTable = GetLookupTableForArray( "Elevation", 1, RGBPoints=[0.2113402634859085, 0.0, 0.0, 1.0, 0.7886597514152527, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039215686275, 0.498039215686275, 0.498039215686275], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_Elevation_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0329900495707989, 0.0, 0.5, 0.0, 0.32995018362999, 1.0, 0.5, 0.0] )


a1_Elevation_PVLookupTable.ScalarOpacityFunction = a1_Elevation_PiecewiseFunction

MeshPartitionFilter2 = MeshPartitionFilter()

DataRepresentation6 = Show()
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation6.SelectionPointFieldDataArrayName = 'Elevation'
DataRepresentation6.ColorArrayName = ('POINT_DATA', 'Elevation')
DataRepresentation6.OpacityArray = 'Elevation'
DataRepresentation6.LookupTable = a1_Elevation_PVLookupTable
DataRepresentation6.ScaleFactor = 0.1
DataRepresentation6.Representation = 'Depth Sort Polygons'
DataRepresentation6.EnablePiston = 1

RenderView1 = GetRenderView()
RenderView1.CameraPosition = [2.386519687194016, 1.1047581669597257, 2.0642281947323715]
RenderView1.CameraParallelScale = 0.8652845525187569
RenderView1.CameraClippingRange = [1.6576813276576574, 5.4729910393285035]


#
# Render N frames for test purposes
#
timingfile=open("timing.txt", "w+")
for num in range(0, 10):
  start =  time()
  RenderView1.ViewTime = num
  Render()
  finish = time()
  timestring = "Render %08f\n" % (finish-start)
  print timestring
  timingfile.write(timestring)
  str = ("snapshot.%04d.png" % num)
  WriteImage(str)

