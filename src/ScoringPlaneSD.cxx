/**
 * @file ScoringPlaneSD.cxx
 * @brief Class defining a basic sensitive detector for scoring planes.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "SimApplication/ScoringPlaneSD.h"

namespace ldmx {

    ScoringPlaneSD::ScoringPlaneSD(G4String name, G4String colName, int subDetID, DetectorID* detID) :
            G4VSensitiveDetector(name),
            detID_(detID),
            subDetID_(subDetID) { 

        // Add the collection name to vector of names.
        this->collectionName.push_back(colName);

        // Register this SD with the manager.
        G4SDManager::GetSDMpointer()->AddNewDetector(this);

        // Set the subdet ID as it will always be the same for every hit.
        detID_->setFieldValue("subdet", subDetID_);
    }

    ScoringPlaneSD::~ScoringPlaneSD() {
        delete detID_;
    }

    G4bool ScoringPlaneSD::ProcessHits(G4Step* step, G4TouchableHistory* history) {

        // Determine if current particle of this step is a Geantino.
        G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
        bool isGeantino = false;
        if (pdef == G4Geantino::Definition() || pdef == G4ChargedGeantino::Definition()) {
            isGeantino = true;
        }

        // Get the edep from the step.
        G4double edep = step->GetTotalEnergyDeposit();

        // Skip steps with no energy dep which come from non-Geantino particles.
        if (edep == 0.0 && !isGeantino) {
            if (verboseLevel > 2) {
                std::cout << "ScoringPlaneSD skipping step with zero edep" << std::endl << std::endl;
            }
            return false;
        }

        // Create a new hit object.
        G4TrackerHit* hit = new G4TrackerHit();

        // Assign track ID for finding the SimParticle in post event processing.
        hit->setTrackID(step->GetTrack()->GetTrackID());

        // Set the edep.
        hit->setEdep(edep);

        // Set the start position.
        G4StepPoint* prePoint = step->GetPreStepPoint();

        // Set the end position.
        G4StepPoint* postPoint = step->GetPostStepPoint();

        G4ThreeVector start = prePoint->GetPosition();
        G4ThreeVector end = postPoint->GetPosition();

        // Set the mid position.
        G4ThreeVector mid = 0.5 * (start + end);
        hit->setPosition(mid.x(), mid.y(), mid.z());

        // Compute path length.
        G4double pathLength = sqrt(pow(start.x() - end.x(), 2) 
                                + pow(start.y() - end.y(), 2) 
                                + pow(start.z() - end.z(), 2));
        hit->setPathLength(pathLength);

        // Set the global time.
        hit->setTime(step->GetTrack()->GetGlobalTime());

        // Set the momentum
        G4ThreeVector p = postPoint->GetMomentum();
        hit->setMomentum(p.x(), p.y(), p.z());

        
        /*
         * Set the 32-bit ID on the hit.
         */
        int cpNumber = prePoint->GetTouchableHandle()->GetCopyNumber();
        std::cout << "Copy number: " << cpNumber << std::endl;
        detID_->setFieldValue(1, cpNumber);
        std::cout << "field value set " << std::endl;
        hit->setID(detID_->pack());
        std::cout << "id packed " << std::endl;
        hit->setLayerID(cpNumber);
        std::cout << "layer id set " << std::endl;

        /*
         * Debug print.
         */
        if (this->verboseLevel > 2) {
            hit->Print();
            std::cout << std::endl;
        }

        // Insert hit into current hits collection.
        hitsCollection_->insert(hit);

        return true;
    }

    void ScoringPlaneSD::Initialize(G4HCofThisEvent* hce) {

        // Setup hits collection and the HC ID.
        hitsCollection_ = new G4TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
        int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        hce->AddHitsCollection(hcID, hitsCollection_);
    }

    void ScoringPlaneSD::EndOfEvent(G4HCofThisEvent*) {

        // Print number of hits.
        if (this->verboseLevel > 0) {
            std::cout << GetName() << " had " << hitsCollection_->entries() << " hits in event" << std::endl;
        }

        // Print each hit in hits collection.
        if (this->verboseLevel > 1) {
            for (unsigned iHit = 0; iHit < hitsCollection_->GetSize(); iHit++) {
                (*hitsCollection_)[iHit]->Print();
            }
        }
    }
}
