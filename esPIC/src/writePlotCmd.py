#!/usr/bin/python

import os

import sys
import getopt
import glob

# ==
# ||
# || Main Program
# ||
# ==

def writePlotCmd(argv):



    # List of all ptcl files
    
    ptclList = glob.glob("ptcl*.plt")

    # Parse the filenames counting the number of PEs for this run

    PElist = []
    timeList = []
    
    for i in range(0,len(ptclList)):
        name = ptclList[i];
        tmp = name.split('_')

        PE  = tmp[1]
        if PE not in PElist:
            PElist.append(PE)

        tmp2 = tmp[2].split(".")
        time = tmp2[0]
        if time not in timeList:
            timeList.append(int(time))


    # Write combined plot command
    
    g = open('pc','w')
#    print >> g , 'reset'
    print >> g , 'set xrange [-.1:1.1]'
    print >> g , 'set yrange [-.1:1.1]'
        
    for time in sorted(timeList):
        print >> g , 'splot \\'
        for PE in sorted(PElist):
            print >> g , '"mesh_' + str(PE) + '_' + str(time) + '.plt" w l , "ptcl_'+str(PE)+'_' + str(time) + '.plt" w d lw 2, \\'
        print >> g , ''
        print >> g , 'pause .01'
    g.close()

    # Write particles only
    
    g = open('pc_ptcl','w')
#    print >> g , 'reset'
    print >> g , 'set xrange [-.1:1.1]'
    print >> g , 'set yrange [-.1:1.1]'
        
    for time in sorted(timeList):
        print >> g , 'plot \\'
        for PE in sorted(PElist):
            print >> g , '"ptcl_'+str(PE)+'_' + str(time) + '.plt" w d lw 2, \\'
        print >> g , ''
        print >> g , 'pause .01'
    g.close()

    # Write combined plot command
    
    g = open('pc_phi','w')
#    print >> g , 'reset'
    print >> g , 'set xrange [-.1:1.1]'
    print >> g , 'set yrange [-.1:1.1]'
        
    for time in sorted(timeList):
        print >> g , 'splot \\'
        for PE in sorted(PElist):
            print >> g , '"phi_' + str(PE) + '_' + str(time) + '.plt" w l , "ptcl_'+str(PE)+'_' + str(time) + '.plt" w d lw 2, \\'
        print >> g , ''
        print >> g , 'pause .01'
    g.close()

    # Write mesh only command
    
    g = open('pc_mesh','w')
#    print >> g , 'reset'
    print >> g , 'set xrange [-.1:1.1]'
    print >> g , 'set yrange [-.1:1.1]'
        
    for time in sorted(timeList):
        print >> g , 'splot \\'
        for PE in sorted(PElist):
            print >> g , '"mesh_' + str(PE) + '_' + str(time) + '.plt" w l , \\'
        print >> g , ''
        print >> g , 'pause .01'
    g.close()


if __name__ == "__main__":
    writePlotCmd(sys.argv[1:])
