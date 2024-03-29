"""Configuration module for dark brem simulation"""

class DarkBremModel() :
    """Storage for parameters of a dark brem model

    All other models should inherit from this class
    in order to keep the correct internal parameters.

    Parameters
    ----------
    name : str
        Name of this dark brem model
    """

    def __init__(self,name) :
        self.name = name

    def __str__(self) :
        string = '%s {'%self.name
        for key in self.__dict__ :
            if key is not self.name :
                string += ' %s=%s'%( key , self.__dict__[key] )
        string += ' }'
        return string

class G4DarkBreM(DarkBremModel) :
    """Configuration for the G4DarkBreM dark brem model

    Parameters
    ----------
    library_path : str
        Path to directory of LHE files containing dark brem events

    Attributes
    ----------
    method : str
        Interpretation method for LHE files
    threshold : float
        Minimum energy [GeV] that electron should have for dark brem to have nonzero xsec
    epsilon : float
        Epsilon for dark brem xsec calculation
    """

    def __init__(self, library_path) :
        super().__init__('vertex_library')
        self.library_path = library_path
        self.method       = 'forward_only'
        self.threshold    = 2.0 #GeV
        self.epsilon      = 0.01

class DMG4(DarkBremModel) :
    """Configuration for the DMG4 model

    Attributes
    ----------
    threshold : float
        Minimum energy [GeV] that electron should have for dark brem to have nonzero xsec
    epsilon : float
        Epsilon for dark brem xsec calculation
    """

    def __init__(self) :
        super().__init__('dmg4')
        self.epsilon      = 0.01
        self.tungsten_target()

    def tungsten_target(self) :
        self.targetA = 183.84
        self.targetZ = 74.0
        self.density = 19.30

    def lead_target(self) :
        self.targetA = 207.2
        self.targetZ = 82.0
        self.density = 11.35

    def copper_target(self) :
        self.targetA = 63.546
        self.targetZ = 29.0
        self.density = 8.960

class DarkBrem:
    """Storage for parameters of dark brem process

    Attributes
    ----------
    ap_mass : float
        Mass of A' in MeV
    global_bias : float
        Biasing factor to apply to the dark brem cross section everywhere (Default: 1.)
    enable : bool
        Should we use the custom Geant4 dark brem process? (Default: No)
    muons : bool
        Do dark brem off muons rather than electrons (Default: False)
    only_one_per_event : bool
        Should we deactivate the process after one dark brem or allow for more than one? (Default: No)
    cache_xsec : bool
        Should we cache the xsec's computed from the model? (Default: yes)
    model : DarkBremModel
        The model that should be use for dark bremming
    """

    def __init__(self) : 
        self.ap_mass            = 0.
        self.global_bias        = 1.
        self.only_one_per_event = False
        self.enable             = False #off by default
        self.muons              = False
        self.cache_xsec         = True
        self.model              = DarkBremModel('UNDEFINED')

    def activate(self, ap_mass, model = None, muons = False) :
        """Activate the dark brem process with the input A' mass [MeV] and dark brem model

        If no dark brem model is given, we do not activate the process
        and only define the A' mass. This allows for some backwards
        compatibility by allowing users to use the LHEPrimaryGenerator
        with A' particles.
        """

        self.ap_mass = ap_mass

        if model is not None :
            if not isinstance(model,DarkBremModel) :
                raise Exception('Dark brem process needs to be configured with an associated DarkBremModel.')
    
            self.enable = True
            self.model  = model

        self.muons = muons

    def __str__(self): 
        """Stringify the DarkBrem configuration

        Returns
        -------
        str
            A human-readable version of all its attributes
        """

        string  = "{ Enabled: %r"%self.enable
        if self.enable :
            string += ", Mass: %.1f MeV"%self.ap_mass
            string += ", Only One Per Event: %r"%self.only_one_per_event
            string += ", Model: %s"%self.model

        string += " }"

        return string
