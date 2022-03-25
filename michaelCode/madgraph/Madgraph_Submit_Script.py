#!/usr/bin/python2

import sys

#Give condor submit file python filename + jobnumber for seeding
if sys.argv[1].isdigit:
    JobNumber = int(sys.argv[1])
if len(sys.argv) < 1:
    sys.exit("Not enough arguments")

offset = 1
MA = 10.0
Mmu = 0.1
#energies = [str(i) for i in range(10,1000,5)]
energies = ["10","20","30","40","50","60","70","80","90","100","120","140","160","180","200","220","240","260","280","300","350","400","450","500","550","600","700","800","900","1000","1100","1200","1400","1600","1800","2000","3000","4000","5000"]

from subprocess import call
#Call order is script, A' mass, beam energy, seed, output folder name, number of events
call(["/local/cms/user/revering/madgraph/MuRunUndecayed.sh", str(MA) ,energies[JobNumber%len(energies)], str(JobNumber+offset), "map_10p0", "10000"])
