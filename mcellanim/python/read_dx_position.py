#!/usr/bin/python

import struct
import sys
import Numeric as N
import array

def extract_molecule_names(file_path, objects):
    array_count =  0
    seek = 0
    for object in objects:
        #if object['label'].startswith('class'):
        #    numbers.append(object['number'])
        if not (object['label'] == 'field' or object['label'] == 'group'):
            array_count += 1
        if object['label'] == 'group':
            seek = object['seek_data'] 
        
    # Extract the names
    names_tmp = {}
    file = open(file_path, 'r')
    #print seek
    file.seek(seek)
    l = file.readline()
    while l:
        # print l.split()
        name = l.split()[1].split('"')[1]
        num = l.split()[2].split('"')[1]
        names_tmp[num] = name
        l = file.readline()
    file.close()

    # Give the name to each molecular object
    names = []
    for object in objects:
        if object['label'].startswith('class'):
            num_tmp = str(int(object['number']) + array_count)
            name = names_tmp[num_tmp]
            names.append(name) 
    return names

def extract_molecule_positions(file_path, objects):
    positions = []
    for object in objects:
        if object['label'].startswith('class'):
            seek = object['seek_data']
            file = open(file_path, 'rb')
            file.seek(seek)
            num_pos = int(object['label'].split('items ')[1].split(' lsb')[0])
            binvalues = array.array('f')
            binvalues.read(file, num_pos*3)
            data = N.array(binvalues, typecode=N.Float)
            data = N.reshape(data, (num_pos,3))
            data =  N.transpose(data)
            #print data
            file.close()
            positions.append(data)
    return positions

def extract_molecules(file_path):
    objects=[]
    file = open(file_path, 'r')
    lines = file.readlines()
    line_count = len(lines)
    file.seek(0)
    for i in range(line_count):
        l = file.readline() 
        if l.startswith('object'):
            object = {}
            # Get the number for each object
            number = l.split()[1].split('"')[1]
            if number.startswith('null'):
                continue
            object['number'] = number
            # Get the label for each object
            label = l.split()[2]
            if label == "class":
                label = 'class' + l.split('class')[1].rstrip()
            object['label'] = label
            # Record the seek to extract the data
            object['seek_data'] = file.tell()
            objects.append(object)
    file.close()
    return objects

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "usage : "+sys.argv[0]+" my_data.dx"
        sys.exit()
    objects = extract_molecules(sys.argv[1])
    #print objects
    positions = extract_molecule_positions(sys.argv[1], objects)
    print positions
    names = extract_molecule_names(sys.argv[1], objects)
    print names
