#include "DMParticleScalar.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleScalar * DMParticleScalar::theInstance = nullptr;

DMParticleScalar* DMParticleScalar::Definition(G4double MassIn, G4double epsilIn)
{
  if( theInstance ) {
    return theInstance;
  }
  const G4String name = "DMParticleScalar";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
                /* Name ..................... */ name,
                /* Mass ..................... */ MassIn,
                /* Decay width .............. */ 0.,
                /* Charge ................... */ 0.*eplus,
                /* 2*spin ................... */ 0,
                /* parity ................... */ 0,
                /* C-conjugation ............ */ 0,
                /* 2*Isospin ................ */ 0,
                /* 2*Isospin3 ............... */ 0,
                /* G-parity ................. */ 0,
                /* type ..................... */ "boson",
                /* lepton number ............ */ 0,
                /* baryon number ............ */ 0,
                /* PDG encoding ............. */ 5400022, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
                /* stable ................... */ true,
                /* lifetime.................. */ 0,
                /* decay table .............. */ NULL,
                /* shortlived ............... */ false,
                /* subType .................. */ "DMParticleScalar",
                /* anti particle encoding ... */ 5400022
            );
  }
  theInstance = reinterpret_cast<DMParticleScalar*>(anInstance);
  return theInstance;
}


DMParticleScalar* DMParticleScalar::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  G4cout << "DMParticleScalar::Definition() : Scalar is not yet defined by the proper method with arguments, exiting" << G4endl;
  exit(1);
}
