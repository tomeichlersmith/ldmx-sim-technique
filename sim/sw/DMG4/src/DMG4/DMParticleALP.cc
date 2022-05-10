#include "DMParticleALP.hh"

#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleALP * DMParticleALP::theInstance = nullptr;

DMParticleALP* DMParticleALP::Definition(G4double MassIn, G4double epsilIn)
{
  if( theInstance ) {
    return theInstance;
  }
  const G4String name = "DMParticleALP";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  G4double WidthIn = 1./(64.*pi)*MassIn*MassIn*MassIn*epsilIn*epsilIn;
  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
                /* Name ..................... */ name,
                /* Mass ..................... */ MassIn,
                /* Decay width .............. */ WidthIn,
                /* Charge ................... */ 0.,
                /* 2*spin ................... */ 0,
                /* parity ................... */ -1,
                /* C-conjugation ............ */ 0,
                /* 2*Isospin ................ */ 0,
                /* 2*Isospin3 ............... */ 0,
                /* G-parity ................. */ 0,
                /* type ..................... */ "boson",
                /* lepton number ............ */ 0,
                /* baryon number ............ */ 0,
                /* PDG encoding ............. */ 5300122, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
                /* stable ................... */ false,
                /* lifetime.................. */ 0,
                /* decay table .............. */ NULL,
                /* shortlived ............... */ false,
                /* subType .................. */ "DMParticleALP",
                /* anti particle encoding ... */ 5300122
            );

    // Life time is given from width
    ((DMParticle*)anInstance)->CalculateLifeTime();

    //create Decay Table
    G4DecayTable* table = new G4DecayTable();

    // create a decay channel
    G4VDecayChannel* mode;
    // ALP -> gamma + gamma
    mode = new G4PhaseSpaceDecayChannel("DMParticleALP", 1., 2, "gamma", "gamma");
    table->Insert(mode);

    anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<DMParticleALP*>(anInstance);
  return theInstance;
}


DMParticleALP* DMParticleALP::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  G4cout << "DMParticleALP::Definition() : ALP is not yet defined by the proper method with arguments, exiting" << G4endl;
  exit(1);
}
