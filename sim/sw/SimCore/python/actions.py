""" actions

actions needed for validation of G4DarkBreM
"""

from LDMX.SimCore import simcfg

class EcalDarkBremFilter(simcfg.UserAction):
    """ Configuration for filtering A' events

    Parameters
    ----------
    minApEnergy : float
        Minimum A' energy to keep the event [MeV]
    """

    def __init__(self,minApEnergy):
        super().__init__('ecal_db_filter','biasing::EcalDarkBremFilter')
        from LDMX.Framework.ldmxcfg import Process
        Process.addLibrary('@CMAKE_INSTALL_PREFIX@/lib/libSimCore_Actions.so')
        self.threshold = minApEnergy

class StepPrinter(simcfg.UserAction) :
    """Print each step of the input track ID

    The default track ID is 1 (the primary particle).

    Parameters
    ----------
    track_id : int, optional
        Geant4 track ID to print each step of
    """

    def __init__(self,track_id=1) :
        super().__init__('print_steps_%s'%track_id,'biasing::utility::StepPrinter')
        from LDMX.Framework.ldmxcfg import Process
        Process.addLibrary('@CMAKE_INSTALL_PREFIX@/lib/libSimCore_Actions.so')
        self.track_id = track_id

class PartialEnergySorter(simcfg.UserAction) :
    """Process particles such that all particles above
    the input threshold are processed first.

    Parameters
    ----------
    thresh : float
        Minimum energy [MeV] to process a track first
    """

    def __init__(self,thresh) :
        super().__init__('sort_above_%dMeV'%thresh,'biasing::utility::PartialEnergySorter')
        from LDMX.Framework.ldmxcfg import Process
        Process.addLibrary('@CMAKE_INSTALL_PREFIX@/lib/libSimCore_Actions.so')
        self.threshold = thresh

