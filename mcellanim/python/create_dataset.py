#!/usr/bin/python

import sys
import os
import numpy
import h5py

from Numeric import *

import read_dx_position

#---------------- Helpfull functions ------------------------------

ERROR_STR= """Error removing %(path)s, %(error)s """

def rmgeneric(path, __func__):
    try:
        __func__(path)
    except OSError, (errno, strerror):
        print ERROR_STR % {'path' : path, 'error': strerror }

def removeall(path):
    if not os.path.isdir(path):
        return
    files=os.listdir(path)
    for x in files:
        fullpath=os.path.join(path, x)
        if os.path.isfile(fullpath):
            f=os.remove
            rmgeneric(fullpath, f)
        elif os.path.isdir(fullpath):
            removeall(fullpath)
            f=os.rmdir
            rmgeneric(fullpath, f)

def create_dir(dir_path):
    if os.path.exists(dir_path):
        removeall(dir_path)
        os.rmdir(dir_path)
    os.mkdir(dir_path)

def add_padding(path, padding):
    files=os.listdir(path)
    for x in files:
        fullpath=os.path.join(path, x)
        if os.path.isfile(fullpath):
            num1 = fullpath.split('.')[-2]
            if len(num1) >= padding:
                continue
            num2 = num1.zfill(padding)
            new_filepath = fullpath.replace('.'+num1+'.', '.'+num2+'.')
            os.system("mv "+fullpath+" "+new_filepath)

#---------------- Main ------------------------------

# Get the input directory with the dx files for each time step
if len(sys.argv) < 3:
    print "usage : ./"+sys.argv[0]+" input/my_data_directory output/out_file.h5part"
    sys.exit()
dx_dir_path = sys.argv[1]
h5part_path = sys.argv[2]

# Create the output directory or overwrite the existing one.
create_dir(h5part_path)

# Add some padding to the file name
add_padding(dx_dir_path, 4)

# Get the list of dx files
#dx_files = sorted([f in os.listdir(dx_dir_path) if os.path.isfile(f)])
dx_files = sorted([dx_dir_path+'/'+f for f in os.listdir(dx_dir_path)])
for f in dx_files :
    if not os.path.isfile(f):
        dx_files.remove(f)
#print dx_files 

# Create the data structure for all the variables (types, positions)
molecules = {}
for path in dx_files:
    print(path)
    objects   = read_dx_position.extract_molecules(path)
    positions = read_dx_position.extract_molecule_positions(path, objects)
    names = read_dx_position.extract_molecule_names(path, objects)
    for i in range(len(names)):
        name = names[i]
        frame_position = positions[i]
        if not molecules.has_key(name):
            molecules[name] = []
        molecules[name].append(frame_position)

#print molecules

# Write the HDF5
print("Writing the output file "+h5part_path)
for molecule in molecules:
    path = h5part_path + "/" + molecule + ".h5part"
    frames = molecules[molecule]
    print path
    f=h5py.File(path,'w')
    for frame_id in range(len(frames)):
        step_grp = f.create_group("Step#"+str(frame_id))
        step_grp.attrs["TimeValue"]=[float(frame_id)]
        for i in range(3):
            label = 'Coords_' + str(i) 
            step_grp.create_dataset(label,data=[float(val) for val in frames[frame_id][i]])
        step_grp.create_dataset('Type',data=zeros(7, Float)) 
