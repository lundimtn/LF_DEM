#!/Applications/VPython-Py2.7/VIDLE.app/Contents/MacOS/Python
# coding=utf-8

# if using Mac OSX, interpreter is
#/Applications/VPython-Py2.7/VIDLE.app/Contents/MacOS/Python

# if using Linux, interpreter is
# /usr/bin/python

import visual
from visual import graph


import sys
import math
from pyLF_DEM_posfile_reading import *
from string import *


def init ():
    global part_nb, phi
    
    if len(sys.argv) == 4 :  # input file 
        input_stream=open(sys.argv[1], 'r')
        part_nb=int(sys.argv[2])
        phi=float(sys.argv[3])

    if len(sys.argv) == 3 : # stdin
        input_stream=sys.stdin
        part_nb=int(sys.argv[1])
        phi=float(sys.argv[2])

    if len(sys.argv) !=3 and len(sys.argv) !=4:
        print ""
        print "use as follows:"
        print "(1) ",sys. argv[0], " N phi < some_input_config"
        print "or" 
        print "(2) ", sys.argv[0], " some_input_config N phi "
        print ""
        sys.exit(1)
    return input_stream
    

# end of init()

def init_visualization():
    global particles, box, scene
    
    half_cell=0.5*pos_stream.lx
#    box=visual.box(pos=(0., 0., 0.), length=pos_stream.cell_lin_size(), width=pos_stream.cell_lin_size(), height=pos_stream.cell_lin_size(), opacity=0.1)
    particles=[ visual.sphere(pos=(half_cell,half_cell,half_cell), radius=0., color=(0.9,0.9,0.9), opacity=1.) for i in range(part_nb) ]
    
def update_visualization():
    global particles
    count=0

    for j in pos_stream.range():

        particles[j].pos=(pos_stream.positions[j,0], pos_stream.positions[j,2], pos_stream.positions[j,1])

        particles[j].color=(0.2, 0.2, 0.5)
        particles[j].radius=pos_stream.radius[j]
        particles[j].material=visual.materials.wood
        particles[j].visible=True

    
# main
part_nb=int()
phi=float()

particles=[]

stream=init()


print ""
print " Plotting visualization for ", stream
print ""

pos_stream=Pos_Stream(stream)
init_visualization()
while pos_stream.get_snapshot():
    update_visualization()
    visual.rate(1000)

