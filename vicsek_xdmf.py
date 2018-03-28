#!/usr/bin/python
# Python script to generate xdmf file.
# Requires configure file >= v2.0
#import argparse
import sys
import getopt
from xdmf import *

##### MAIN ####
#parser = argparse.ArgumentParser(description='xdmf converter')
#parser.add_argument('-f', '--file', action="store", help='file name',
#                    dest="file")
#parser.add_argument('-t', '--time', action="store", help='time frames',
#                    dest="time", type=int)
#parser.add_argument('-n', '--num',  action="store", help='number particles',
#                    dest="nump")
#results    = parser.parse_args()
try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:f:t:n:', ["dir", "file", "time", "num"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

h5dir = h5file = max_frames = nump = None
for o, a in opts:
    if o in ("-d", "--dir"):
        h5dir      = str(a)
    elif o in ("-f", "--file"):
        h5file     = str(a)
    elif o in ("-t", "--time"):
        max_frames = int(a)
    elif o in ("-n", "--num"):
        nump       = int(a)
    else:
        assert False, "unknown option"
        
if h5dir == None or h5file == None or max_frames == None or nump == None:
    parser.print_help()
    sys.exit(2)

tag  = h5file.replace('.h5','')
time = [0, 1.0, max_frames]
px   = init(time, "Vicsek Particle Trajectory")
for i in range(max_frames):
    loc = h5file+':/trj/t_'+str(i)+'/'
    conf= newPoly2D(px, nump, loc+'pos', 'pos')
    addScalar(px, conf, loc+'num', 'num', dtype=int4)
    addScalar(px, conf, loc+'num_v', 'num_v', dtype=int4)
    addScalar(px, conf, loc+'num_xi', 'num_xi', dtype=int4)
    addScalar(px, conf, loc+'phi', 'phi')
    addVector2D(px, conf, loc+'vel', 'vel')
    freePoly(px)
free(px, h5dir+"/"+tag+'.xdmf')

bx   = init([0, 1.0, 1], "Vicsek Box Data")
loc  = h5file+':/sys/'
conf = newGrid2D(bx, [2, 2], [loc+'box_x', loc+'box_y'], 'box')
freePoly(bx)
free(bx, h5dir+"/"+tag+'_box.xdmf')
