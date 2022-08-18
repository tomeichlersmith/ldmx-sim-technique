import sys
import os
import argparse

usage = "ldmx fire %s"%(sys.argv[0])
parser = argparse.ArgumentParser(usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Both method options
parser.add_argument('--out_dir',default=os.getcwd(),help='Directory to output data files')
parser.add_argument("-p", "--pause",dest="pause",default=False,action='store_true',
        help='Print the process and pause to continue processing.')
parser.add_argument("-v","--verbose",dest="verbose",default=False,action='store_true',
        help='Print periodic progress messages.')
parser.add_argument('-d','--depth',default=100.,type=float,
        help='Depth of material hunk [mm] to simulate inside of.')
parser.add_argument('--max_events',default=50000,type=int,
        help='Maximum number of events to produce')

dbmethod_parser = parser.add_subparsers(title='method',
    help='Choose a method for simulating dark brem')

dmg4_parser = dbmethod_parser.add_parser('dmg4',help='Use DMG4 to calculate xsec and kinematics')
dmg4_parser.set_defaults(method = 'dmg4')
dmg4_parser.add_argument("run_number",default=1,type=int,
        help='Run number controlling rand numbers')
dmg4_parser.add_argument('-m','--material',default='tungsten',
        choices=['tungsten','silicon','brass','lead'],
        help='Material to use as target in simulation.')
dmg4_parser.add_argument('--particle',default='electron',
        choices=['electron','muon'],
        help='Particle to be the primary in the simulation.')
dmg4_parser.add_argument('--primary_energy',default=4.,type=float,
        help='Energy of primary particle [GeV].')
dmg4_parser.add_argument('--ap_mass',default=1000.,type=float,
        help='A prime mass [MeV].')

g4db_parser = dbmethod_parser.add_parser('g4db',help='Use G4DarkBreM to calculate xsec and kinematics')
g4db_parser.set_defaults(method = 'g4db')
g4db_parser.add_argument("db_lib",type=str,default=None,
        help="Archive or directory for a dark brem event library to use for the model.")

arg = parser.parse_args()

hunk_transverse = 500 #mm

from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('db')
p.maxEvents = arg.max_events
p.maxTriesPerEvent = 1000

if arg.method == 'g4db' :
    # Dark Brem Event Library
    if arg.db_lib.endswith('.tar.gz') :
        # 1) Unpack the archive
        import tarfile
        with tarfile.open(arg.db_lib,"r:gz") as ar :
            ar.extractall()
    
        # 2) Define path to library
        #   extracting the library puts the directory in the current working directory
        #   so we just need the basename
        db_event_lib_path = os.path.basename(arg.db_lib).replace('.tar.gz','')
    else :
        # We were given the directory of the event library,
        #   so we don't have to do anything else
        db_event_lib_path = arg.db_lib
    
    # need to remove trailing slash so we can deduce the parameters from the library name
    if db_event_lib_path.endswith('/') :
        db_event_lib_path = db_event_lib_path[:-1]
    
    # Get A' mass and run number from the dark brem library name
    lib_parameters = os.path.basename(db_event_lib_path).split('_')
    particle = lib_parameters[0]
    material = lib_parameters[1]
    ap_mass = float(lib_parameters[lib_parameters.index('mA')+1])*1000. # MeV for A' mass
    run_num = int(lib_parameters[lib_parameters.index('run')+1]) # run number
    max_e = float(lib_parameters[lib_parameters.index('MaxE')+1]) # GeV for ParticleGun
    min_e = float(lib_parameters[lib_parameters.index('MinE')+1])*1000. # MeV for user actions
    # using copper as library for brass
    if material == 'copper' :
        material = 'brass'
elif arg.method == 'dmg4' :
    particle = arg.particle
    material = arg.material
    ap_mass = arg.ap_mass
    run_num = arg.run_number
    max_e = arg.primary_energy
    min_e = 2.*arg.ap_mass

if not os.path.isdir(arg.out_dir) :
    os.makedirs(arg.out_dir)

p.outputFiles = ['/dev/null']

p.histogramFile = f'{arg.out_dir}/ntuple_{arg.method}_{particle}_{material}_depthmm_{arg.depth}_mAMeV_{int(ap_mass)}_events_{p.maxEvents}_run_{run_num}.root'
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
primary = generators.gun(f'primary_{particle}')
if particle == 'electron' :
    primary.particle = 'e-'
elif particle == 'muon' :
    primary.particle = 'mu-'
else :
    raise Exception(f'Unrecognized primary particle {particle}')

primary.energy = max_e
primary.direction = [ 0., 0., 1. ] #unitless
primary.position = [ 0., 0., -1. ] #mm
sim.generators = [ primary ]
sim.time_shift_primaries = False

# use detector.gdml file in current directory
from LDMX.SimCore.write_detector import write
sim.detector = write(p.histogramFile.replace("root","gdml").replace('ntuple','geo'),
    material,arg.depth,hunk_transverse)

#Activiate dark bremming with a certain A' mass and LHE library
from LDMX.G4DarkBreM import dark_brem
if arg.method == 'g4db' :
    db_model = dark_brem.G4DarkBreM(db_event_lib_path)
    db_model.threshold = 2*min_e/1000. #GeV - minimum energy electron needs to have to dark brem
    db_model.epsilon   = 0.01 #decrease epsilon from one to help with Geant4 biasing calculations
else :
    db_model = dark_brem.DMG4()
    db_model.epsilon = 0.0001 # default in DMG4

sim.dark_brem.activate( ap_mass , db_model , muons = (primary.particle == 'mu-'))
sim.dark_brem.only_one_per_event = True

#Biasing dark brem up inside of the ecal volumes
from math import log10
#need a higher power for the higher mass A'
mass_power = max(log10(sim.dark_brem.ap_mass),2.)
if particle == 'muon' :
    mass_power = 1

from LDMX.SimCore import bias_operators
sim.biasing_operators = [ 
        bias_operators.DarkBrem('hunk',True,
          sim.dark_brem.ap_mass**mass_power / db_model.epsilon**2,
          particle = primary.particle)
        ]

from LDMX.Biasing import filters
from LDMX.Biasing import util
sim.actions = [ 
    util.PartialEnergySorter(min_e), 
    filters.EcalDarkBremFilter(0.), 
#    util.StepPrinter(1) 
    ]

p.sequence = [
    sim,
    ldmxcfg.Analyzer('dbint','dqm::NtuplizeDarkBremInteraction','Biasing')
    ]

if arg.pause :
    p.pause()
