#ifndef SIMCORE_DETECTORCONSTRUCTION_H
#define SIMCORE_DETECTORCONSTRUCTION_H

//---< Geant4 >---//
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VUserDetectorConstruction.hh"

//---< Framework >---//
#include "Framework/Configure/Parameters.h"

namespace simcore {

/**
 * @class DetectorConstruction
 * @brief Implements the Geant4 detector construction
 *
 * @note
 * This class reads in a detector description from a GDML file
 * using the basic <i>G4GDMLParser</i> and instantiates supplemental
 * information using the AuxInfoReader.
 *
 * @see AuxInfoReader
 */
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  /**
   * Constructor.
   *
   * @param parser Parser used to parse the geometry into memory.
   * @param parameters The parameters used to configure this class.
   * @param ci The conditions needed to build the detector.
   */
  DetectorConstruction(framework::config::Parameters &parameters);

  /**
   * Class destructor.
   */
  ~DetectorConstruction() = default;

  /**
   * Construct the detector.
   * @return The top volume of the detector.
   */
  G4VPhysicalVolume *Construct();

  /**
   */
  void ConstructSDandField();

 private:
  std::string gdml_file_;
  bool validate_;
};  // DetectorConstruction
}  // namespace simcore

#endif  // SIMCORE_DETECTORCONSTRUCTION_H
