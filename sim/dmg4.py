import sys
import os
import argparse

usage = "ldmx fire %s"%(sys.argv[0])
parser = argparse.ArgumentParser(usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-p", "--pause",dest="pause",default=False,action='store_true',
        help='Print the process and pause to continue processing.')
parser.add_argument("-v","--verbose",dest="verbose",default=False,action='store_true',
        help='Print periodic progress messages.')

parser.add_argument('-d','--depth',default=100.,type=float,
        help='Depth of material hunk [mm] to simulate inside of.')
parser.add_argument('-m','--material',default='tungsten',
        choices=['tungsten','silicon'],
        help='Material to use as target in simulation.')
parser.add_argument('--particle',default='e-',
        help='Particle to be the primary in the simulation.')
parser.add_argument('--primary_energy',default=4.,type=float,
        help='Energy of primary particle [GeV].')

parser.add_argument('--out_dir',default=os.getcwd(),
        help='Directory to output event file.')

arg = parser.parse_args()

hunk_transverse = 500 #mm

from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('db')
p.maxEvents = 1#0
p.maxTriesPerEvent = 10 #00

# Get A' mass and run number from the dark brem library name
ap_mass = 10.
run_num = 1

if not os.path.isdir(arg.out_dir) :
    os.makedirs(arg.out_dir)

p.outputFiles = [
        f'{arg.out_dir}/{arg.particle}_{arg.material}_depthmm_{arg.depth}_mAMeV_{int(ap_mass)}_events_{p.maxEvents}_run_{run_num}.root'
        ]

p.histogramFile = f'{arg.out_dir}/ntuple_{os.path.basename(p.outputFiles[0])}'
p.run = run_num

from LDMX.SimCore import simulator
sim = simulator.simulator( "dark_brem_%sMeV" % str(ap_mass) )
sim.description = "Dark Brem Process Testing and Validation"

if arg.verbose :
    sim.validate_detector = True
    sim.verbosity = 2
    sim.preInitCommands = ['/run/verbose 2']
    p.termLogLevel = 0
    p.logFrequency = 1

from LDMX.SimCore import generators
primary = generators.gun(f'primary_{arg.particle}')
primary.particle = arg.particle
primary.energy = arg.primary_energy
primary.direction = [ 0., 0., 1. ] #unitless
primary.position = [ 0., 0., -1. ] #mm
sim.generators = [ primary ]
sim.time_shift_primaries = False

# use detector.gdml file in current directory
from LDMX.Detectors.write import write
sim.detector = write('.write_detector.gdml',arg.material,hunk_transverse,arg.depth)

#Activiate dark bremming with a certain A' mass and LHE library
from LDMX.SimCore import dark_brem
db_model = dark_brem.DMG4Model()
db_model.threshold = 2. #GeV - minimum energy electron needs to have to dark brem
sim.dark_brem.activate( ap_mass , db_model )

#Biasing dark brem up inside of the ecal volumes
from math import log10
#need a higher power for the higher mass A'
mass_power = max(log10(sim.dark_brem.ap_mass),2.)

from LDMX.SimCore import bias_operators
sim.biasing_operators = [ 
        bias_operators.DarkBrem('hunk',True,1e300) #sim.dark_brem.ap_mass**mass_power / 0.01**2)
        ]

from LDMX.Biasing import filters
from LDMX.Biasing import util
sim.actions = [
        #Make sure all particles above 2GeV are processed first
        util.PartialEnergySorter(2000.),
        #Only keep events when a dark brem happens in the Ecal (the hunk in this case)
        filters.EcalDarkBremFilter(2000.)
]

p.sequence = [sim]

if arg.pause :
    p.pause()
