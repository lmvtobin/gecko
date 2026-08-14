"""
* * * * * * * * * * * * * * * * * * * * *
*                           o---------x *
*     G e c K o         x-o |  #  w w | *
*                     o-|-x |       w | *
*       Grid          x-|---|---o     | *
*       Kirby     x-----o   |   |   # | *
*                 |     x---|---|-----o *
*             x---o     |   |   | x-o   *
*             |     o---|---x o-|-|-x   *
*         o-x |     |   |     | x-o     *
*       x-|-|-o x---|---o     |         *
*       o-x |   |   |     o---x         *
* o---------|---|---x     |             *
* |         |   |   o-----x             *
* |     o---|---|---|-x                 *
* |     |   |   | x-|-o        By       *
* |   x-|---o   | o-x         Lucy      *
* | o-|-x       |             Tobin     *
* x-|-o         |                       *
*   x-----------o                       *
* * * * * * * * * * * * * * * * * * * * *
"""

# Some additional python utilities for reading in sample data files
# In future, this will be replaced by a proper python interface to geckocore

import os.path
from regina import *

def read_sample_data(filename):
    pairs = []
    with open(filename,'r') as f:
        for line in f:
            sigs_string,count,grids_string = line.split(" : ")
            sigs = sigs_string.split(" ")
            grids = []
            for s in grids_string.split(" "):
                if s != "\n":
                    x_string,o_string = s[1:-1].split("][")
                    x = [int(k) for k in x_string.split(",")]
                    o = [int(k) for k in o_string.split(",")]
                    grids.append((x,o))
            pairs.append([sigs,grids])
    return pairs

def read_volume_data(filename):
    hyp = []
    nonhyp = []
    with open(filename,'r') as f:
        reading_hyp = True
        for line in f:
            if line.startswith("0x") or line.startswith("-0x"):
                volstr,gridstr = line.split(" ")
                vol = float.fromhex(volstr)
                xstr,ostr = gridstr.rstrip()[1:-1].split("][")
                x = [int(c) for c in xstr.split(",")]
                o = [int(c) for c in ostr.split(",")]
                if reading_hyp:
                    hyp.append((vol,(x,o)))
                else:
                    nonhyp.append((vol,(x,o)))
            elif "likely non-hyperbolic:" in line:
                reading_hyp = False
    return (hyp,nonhyp)


def read_form_data(filename):
    pairs = []
    with open(filename,'r') as f:
        for line in f:
            form_string,count,grids_string = line.split(" : ")
            form = [int(s) for s in form_string.split(" ")]
            grids = []
            for s in grids_string.split(" "):
                if s != "\n":
                    x_string,o_string,assignment = s.split("]")
                    x_string = x_string[1:]
                    o_string = o_string[1:]
                    x = [int(k) for k in x_string.split(",")]
                    o = [int(k) for k in o_string.split(",")]
                    grids.append((x,o,int(assignment)))
            pairs.append([form,grids])
    return pairs

def is_split(tri):
    
    if tri.countBoundaryComponents()<2:
        return False

    # TODO: easy preconditions which confirm not split (e.g. linking number)

    # link is split iff complement has a vertex normal surface which is a sphere
    # and separates two boundary components
    # I believe cutAlong() should always produce real boundaries for newly created
    # boundary comps, so components should contain one of the original boundaries
    # iff they have an ideal vertex
    # TODO: verify this is actually true
    surfaces = NormalSurfaces(tri,NormalCoords.Standard)
    for s in surfaces:
        if s.eulerChar()==2:
            comps = s.cutAlong().components()
            if len(comps)==2 and comps[0].isIdeal() and comps[1].isIdeal():
                return True

    return False
