#!/usr/bin/python

import sys
import os
import numpy
import h5py

# Get the input directory with the ascii files for each time step
if len(sys.argv) < 3:
    print "usage : ./"+sys.argv[0]+" input/my_data_directory output/out_file.h5part"
    sys.exit()
ascii_dir_path = sys.argv[1]
h5part_path = sys.argv[2]

# Get the list of ascii files
#ascii_files = sorted([f in os.listdir(ascii_dir_path) if os.path.isfile(f)])
ascii_files = sorted([ascii_dir_path+'/'+f for f in os.listdir(ascii_dir_path)])

# Set three name of the variables and the indice in the ascii data files
labels = [['Coords_0',1],['Coords_1',2],['Coords_2',3],['Type',0]]

# Set the separator between the values
separator = ' '

# Check that all the files are consistent
f = open(ascii_files[0])
lines = f.readlines()
particule_count = len(lines)
field_count = len(labels)
#field_count = len(lines[0])
f.close()
for path in ascii_files:
    f = open(path)
    lines = f.readlines()
    if not len(lines) == particule_count:
        print("File "+path+" : particule count = "+str(len(lines))+\
              " instead of "+str(particule_count))
        f.close()
        exit()
    if not len(lines[0].split()) == field_count:
        print("File "+path+" : field count = "+str(len(lines[0].split()))+\
              " instead of "+str(field_count))
        f.close()
        exit()
    f.close()

# Create the data structure for all the variables (types, positions)
frames = []
for path in ascii_files:
    print(path)
    frame = numpy.loadtxt(path, unpack=True)
    frames.append(frame)     
    f.close()

# Write the HDF5
print("Writing the output file "+h5part_path)
f=h5py.File(h5part_path,'w')
#root_grp = f.create_group('/')
for frame_id in range(len(frames)):
    #step_grp = root_grp.create_group("Step#"+str(frame_id))
    step_grp = f.create_group("Step#"+str(frame_id))
    step_grp.attrs["TimeValue"]=[float(frame_id)]
    for label in labels:
        index = label[1]
        label_frame_data = [float(val) for val in frames[frame_id][index,:]]
        step_grp.create_dataset(label[0],data=label_frame_data)
