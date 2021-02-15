#include "SimCore/UserEventInformation.h"

/*~~~~~~~~~~~~~~~~*/
/*   C++ StdLib   */
/*~~~~~~~~~~~~~~~~*/
#include <iostream>

namespace simcore {

UserEventInformation::UserEventInformation() {}

UserEventInformation::~UserEventInformation() {}

void UserEventInformation::Print() const {
  std::cout << "Event weight: " << weight_ << "\n"
            << "Brem candidate count: " << bremCandidateCount_ << "\n"
            << std::endl;
}
}  // namespace simcore
