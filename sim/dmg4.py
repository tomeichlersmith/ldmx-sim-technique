import sys
import os
import argparse

usage = "ldmx fire %s"%(sys.argv[0])
parser = argparse.ArgumentParser(usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("run_number",default=1,type=int,
    help='Run number controlling rand numbers')
parser.add_argument("-p", "--pause",dest="pause",default=False,action='store_true',
        help='Print the process and pause to continue processing.')
parser.add_argument("-v","--verbose",dest="verbose",default=False,action='store_true',
        help='Print periodic progress messages.')

parser.add_argument('-d','--depth',default=100.,type=float,
        help='Depth of material hunk [mm] to simulate inside of.')
parser.add_argument('-m','--material',default='tungsten',
        choices=['tungsten','silicon','brass'],
        help='Material to use as target in simulation.')
parser.add_argument('--particle',default='e-',
        choices=['e-','mu-'],
        help='Particle to be the primary in the simulation.')
parser.add_argument('--primary_energy',default=4.,type=float,
        help='Energy of primary particle [GeV].')
parser.add_argument('--ap_mass',default=1000.,type=float,
        help='A prime mass [MeV].')

parser.add_argument('--out_dir',default=os.getcwd(),
        help='Directory to output event file.')

arg = parser.parse_args()

hunk_transverse = 500 #mm

from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('db')
p.maxEvents = 50000
p.maxTriesPerEvent = 1000

if not os.path.isdir(arg.out_dir) :
    os.makedirs(arg.out_dir)

particle = 'electron'
if arg.particle == 'mu-' :
    particle = 'muon'

p.outputFiles = [
        f'{arg.out_dir}/{particle}_{arg.material}_depthmm_{arg.depth}_mAMeV_{int(arg.ap_mass)}_events_{p.maxEvents}_run_{arg.run_number}.root'
        ]

p.histogramFile = f'{arg.out_dir}/ntuple_{os.path.basename(p.outputFiles[0])}'
p.run = arg.run_number

from LDMX.SimCore import simulator
sim = simulator.simulator( "dark_brem_%sMeV" % str(arg.ap_mass) )
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
sim.detector = write('.write_detector.gdml',arg.material,arg.depth,hunk_transverse)

#Activiate dark bremming with a certain A' mass and LHE library
from LDMX.SimCore import dark_brem
db_model = dark_brem.DMG4Model()
sim.dark_brem.activate( arg.ap_mass , db_model , muons = (primary.particle == 'mu-'))
sim.dark_brem.only_one_per_event = True

#Biasing dark brem up inside of the ecal volumes
from math import log10
#need a higher power for the higher mass A'
mass_power = max(log10(sim.dark_brem.ap_mass),2.)

from LDMX.SimCore import bias_operators
bias = 0.01**2
if primary.particle == 'mu-' :
    bias = 0.01

sim.biasing_operators = [ 
        bias_operators.DarkBrem('hunk',True,bias,particle = primary.particle)
        ]

from LDMX.Biasing import filters
from LDMX.Biasing import util
sim.actions = [
        #Make sure all particles above 2GeV are processed first
        util.PartialEnergySorter(2000.),
        #Only keep events when a dark brem happens in the Ecal (the hunk in this case)
        filters.EcalDarkBremFilter(0.)
]

p.sequence = [
    sim,
    ldmxcfg.Analyzer('dbint','dqm::NtuplizeDarkBremInteraction','Biasing')
    ]

if arg.pause :
    p.pause()
