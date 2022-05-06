#include "DMParticleXScalar.hh"

#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleXScalar * DMParticleXScalar::theInstance = nullptr;

DMParticleXScalar* DMParticleXScalar::Definition(G4double MassIn, G4double epsilIn)
{
  if( theInstance ) {
    return theInstance;
  }
  const G4String name = "DMParticleXScalar";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  G4double RatioEA2 = electron_mass_c2*electron_mass_c2/(MassIn*MassIn);
  G4double WidthIn = (1./2.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.-4.*RatioEA2);
  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
                /* Name ..................... */ name,
                /* Mass ..................... */ MassIn,
                /* Decay width .............. */ WidthIn,
                /* Charge ................... */ 0.*eplus,
                /* 2*spin ................... */ 2,
                /* parity ................... */ +1,
                /* C-conjugation ............ */ 0,
                /* 2*Isospin ................ */ 0,
                /* 2*Isospin3 ............... */ 0,
                /* G-parity ................. */ 0,
                /* type ..................... */ "boson",
                /* lepton number ............ */ 0,
                /* baryon number ............ */ 0,
                /* PDG encoding ............. */ 5400122, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
                /* stable ................... */ false,
                /* lifetime.................. */ 0,
                /* decay table .............. */ NULL,
                /* shortlived ............... */ false,
                /* subType .................. */ "DMParticleXScalar",
                /* anti particle encoding ... */ 5400122
            );

    // Life time is given from width
    ((DMParticle*)anInstance)->CalculateLifeTime();

    //create Decay Table
    G4DecayTable* table = new G4DecayTable();

    // create a decay channel
    G4VDecayChannel* mode;
    // X -> e+ + e-
    mode = new G4PhaseSpaceDecayChannel("DMParticleXScalar", 1., 2, "e-", "e+");
    table->Insert(mode);

    anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<DMParticleXScalar*>(anInstance);
  return theInstance;
}


DMParticleXScalar* DMParticleXScalar::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  G4cout << "DMParticleXScalar::Definition() : XScalar is not yet defined by the proper method with arguments, exiting" << G4endl;
  exit(1);
}
