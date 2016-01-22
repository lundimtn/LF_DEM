#!/usr/bin/env python
#
# yapgen.py, part of LF_DEM
#
# Purpose: yapgen.py is a script to generate a Yaplot visualization
#          of a LF_DEM simulation from the par_* and int_* files generated
#          by LF_DEM.
#
#  Romain Mari, 2015

import sys
import pandas as pd
import numpy as np
import lfdem_file as lf


def read_data(posfile, intfile):
    """
    Purpose:
        Read a par_ and int_ file.

    Returning values:
        pos_frames: a list of snapshots from the par_ file
        int_frames: a list of snapshots from the int_ file
        strains: the associated strains
        shear_rates: the associated strain rates
    """
    field_nb = 15

    # load with pandas read_table
    pos_frames, strains, shear_rates = lf.read_snapshot_file(posfile, field_nb)

    # load with pandas read_table
    field_nb = 17
    int_frames, strains, shear_rates = lf.read_snapshot_file(intfile, field_nb)
#    print pos_frames[:3], int_frames[:3]
    return pos_frames, int_frames, strains, shear_rates

def hfill_array(cmd_array):
    """
        Purpose:
            From an array of strings with shape (N,m), with m<=7,
            get an array of strings with shape (N,7).
    """
    fill_nb = 7 - cmd_array.shape[1]
    height = cmd_array.shape[0]
    filler = np.tile([''], (height,fill_nb))
    return np.column_stack((cmd_array,filler))

def cmd(switch_type, switch_value):
    """
        Purpose:
            Get an array of strings for commands of type switch_type taking values
            switch_value.
    """
    sval = np.array(switch_value, dtype = np.str,ndmin=1)
    switch_cmd = np.empty(sval.shape[0], dtype = np.str)
    switch_cmd[:] = switch_type
    cmd = np.column_stack((switch_cmd, sval))
    return hfill_array(cmd)

def add_cmd(yap_array, switch_type, switch_value):
    """
        Purpose:
            Append to yap_array an array of strings for commands of type switch_type taking values
            switch_value.
    """
    return np.row_stack((yap_array, cmd(switch_type, switch_value)))

def layer_switch(value):
    return cmd('y', value)
def add_layer_switch(yap_array, value):
    return add_cmd(yap_array, 'y', value)

def color_switch(value):
    return cmd('@', value)
def add_color_switch(yap_array, value):
    return add_cmd(yap_array, '@', value)

def radius_switch(value):
    return cmd('r', value)
def add_radius_switch(yap_array, value):
    return add_cmd(yap_array, 'r', value)


def pair_cmd_and_switch(cmd, switch):
    """
    Purpose:
        Common use case: you want to change state (e.g. width) for every object.
        You can do that in an array-like fashion, generating cmd and switch arrays separately,
        and blending them afterwards. This is what this function is for.
    """
    return np.reshape(np.column_stack((switch,cmd)),(2*switch.shape[0], switch.shape[1]))


def get_particles_yaparray(pos, rad):
    """
        Get yaplot commands (as an aray of strinfs) to display circles
        for each particle defined in p.
        p must contain positions and radii of particles.
    """

    yap_out = layer_switch(3)
    yap_out = add_color_switch(yap_out,3)

    particle_circle_positions = cmd('c', pos)
    particle_circle_radius = cmd('r', rad)
    particle_out = pair_cmd_and_switch(particle_circle_positions, particle_circle_radius)
    yap_out = np.row_stack((yap_out, particle_out))

    return yap_out

def get_interactions_yaparray(r1r2, thicknesses):
    """
        Get yaplot commands (as an aray of strings) to display sticks
        for each interactions with end points in r1r2.
        thicknesses contains the width of the sticks.
    """
    yap_out = layer_switch(2)
    yap_out = add_color_switch(yap_out,4)

    interaction_sticks = cmd('s', r1r2.astype(np.str))
    interaction_widths = cmd('r', thicknesses)
    interaction_out = pair_cmd_and_switch(interaction_sticks, interaction_widths)
    yap_out = np.row_stack((yap_out, interaction_out))
    return yap_out

def get_interaction_end_points(f,p):
    """
        For each interaction in f, get the position of the particles involved.
        Positions of every particle given in p.

        Returns an array containing x1,y1,z1,x2,y2,z2 for each interaction.
    """
    # for each interaction: the particle indices
    part1 = f[:,0].astype(np.int)
    part2 = f[:,1].astype(np.int)

    # for each interaction: the particle positions
    r1 = p[part1,2:5].astype(np.float)
    r2 = p[part2,2:5].astype(np.float)

    return np.hstack((r1, r2))

def filter_interactions_crossing_PBC(f,r1r2):
    """
        Exclude interactions across the boundaries
    """
    r1 = r1r2[:,:3]
    r2 = r1r2[:,3:]
    keep = np.linalg.norm(r2-r1,axis=1) < 4.
    r1r2 = r1r2[keep]
    f = f[keep]
    return f, r1r2



def snaps2yap(pos_fname, force_factor):
    forces_fname = pos_fname.replace("par_", "int_")
    positions, forces, strains, shear_rates = read_data(pos_fname, forces_fname)

    yap_filename = pos_fname.replace("par_", "y_")
    yap_file = open(yap_filename,'wb')

    nb_of_frames = len(strains)
    i = 0
    for f,p,strain,rate in zip(forces,positions, strains, shear_rates):
        r1r2 = get_interaction_end_points(f,p)
        f, r1r2 = filter_interactions_crossing_PBC(f,r1r2)

        # display a line joining the center of interacting particles
        # with a thickness proportional to the normal force
        normal_forces = (f[:,7]+f[:,11]+f[:,15]).astype(np.float)
        normal_forces = force_factor*np.abs(normal_forces) # to convert the force to a thickness. case-by-case.
        yap_out = get_interactions_yaparray(r1r2, normal_forces)

        # display a circle for every particle
        pos = p[:,2:5].astype(np.float)
        rad = p[:,1].astype(np.float)
        yap_out = np.row_stack((yap_out, get_particles_yaparray(pos, rad)))

        yap_out = add_layer_switch(yap_out, 5)

        # display strain
        yap_out = np.row_stack((yap_out, ['t',str(10.),str(0.),str(10.),'strain='+str(strain),'','']))

        # output
        np.savetxt(yap_file, yap_out, fmt="%s "*7)
        yap_file.write("\n".encode('utf-8'))
        i += 1
        print("\rframe "+str(i)+"/"+str(nb_of_frames),end="",flush=True)
    yap_file.close()

def conf2yap(conf_fname):
    yap_filename = pos_fname.replace(".dat", ".yap")
    yap_file = open(yap_filename,'wb')

    positions, radii, meta = lf.read_conf_file(conf_fname)
    positions[:,0] -= float(meta['lx'])/2
    positions[:,1] -= float(meta['ly'])/2
    positions[:,2] -= float(meta['lz'])/2

    yap_out = get_particles_yaparray(positions, radii)
    # print(yap_out)
    np.savetxt(yap_file, yap_out, fmt="%s "*7)
    yap_file.write("\n".encode('utf-8'))
    yap_file.close()

if len(sys.argv) < 2:
    print(sys.argv[0], " par_or_conf_file [force_factor]\n")
    exit(1)

pos_fname = sys.argv[1]

if pos_fname.find("par_") > -1:
    if len(sys.argv)>2:
        force_factor = float(sys.argv[2])
    else:
        force_factor = 0.01
    snaps2yap(pos_fname, force_factor)
else:
    conf2yap(pos_fname)
