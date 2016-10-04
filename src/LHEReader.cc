#include "SimApplication/LHEReader.h"

// STL
#include <iostream>
#include <stdexcept>

LHEReader::LHEReader(std::string& filename) {
    std::cout << "Opening LHE file " << filename << std::endl;
    ifs.open(filename.c_str(), std::ifstream::in);
}

LHEReader::~LHEReader() {
    ifs.close();
}

LHEEvent* LHEReader::readNextEvent() {

    std::string line;
    bool foundEventElement = false;
    while (getline(ifs, line)) {
        if (line == "<event>") {
            foundEventElement = true;
            break;
        }
    }

    if (!foundEventElement) {
        std::cerr << "WARNING: No next <event> element was found by the LHE reader." << std::endl;
        return NULL;
    }

    getline(ifs, line);

    LHEEvent* nextEvent = new LHEEvent(line);

    std::cout << "  NUP: " << nextEvent->getNUP() << ", IDPRUP: " << nextEvent->getIDPRUP()
            << ", XWGTUP: " << nextEvent->getXWGTUP() << ", SCALUP: " << nextEvent->getSCALUP()
            << ", AQEDUP: " << nextEvent->getAQEDUP() << ", AQCDUP: " << nextEvent->getAQCDUP()
            << std::endl;

    while (getline(ifs, line)) {
        if (line == "</event>") {
            std::cout << "Found end event element." << std::endl;
            break;
        }
        std::cout << "Read LHE record: " << line << std::endl;
        LHEParticle* particle = new LHEParticle(line);

        std::cout << " "
                << "IDUP: " << particle->getIDUP()
                << ", ISTUP: " << particle->getISTUP()
                << ", MOTHUP[0]: " << particle->getMOTHUP(0)
                << ", MOTHUP[1]: " << particle->getMOTHUP(1)
                << ", ICOLUP[0]: " << particle->getICOLUP(0)
                << ", ICOLUP[1]: " << particle->getICOLUP(1)
                << ", PUP[0]: " << particle->getPUP(0)
                << ", PUP[1]: " << particle->getPUP(1)
                << ", PUP[2]: " << particle->getPUP(2)
                << ", PUP[3]: " << particle->getPUP(3)
                << ", PUP[4]: " << particle->getPUP(4)
                << ", VTIMUP: " << particle->getVTIMUP()
                << ", SPINUP: " << particle->getSPINUP()
                << std::endl;

        nextEvent->addParticle(particle);
    }

    const std::vector<LHEParticle*>& particles = nextEvent->getParticles();
    int particleIndex = 0;
    for (std::vector<LHEParticle*>::const_iterator it = particles.begin();
            it != particles.end(); it++) {
        LHEParticle* particle = (*it);
        if (particle->getMOTHUP(0) != 0) {
            int mother1 = particle->getMOTHUP(0);
            int mother2 = particle->getMOTHUP(1);
            if (mother1 > 0) {
                std::cout << "  Assigning mother particle " << mother1 << " to particle at index " << particleIndex << std::endl;
                particle->setMother(0, particles[mother1 - 1]);
            }
            if (mother2 > 0) {
                std::cout << "  Assigning mother particle " << mother2 << " to particle at index " << particleIndex << std::endl;
                particle->setMother(1, particles[mother2 - 1]);
            }
        }
        ++particleIndex;
    }

    return nextEvent;
}
