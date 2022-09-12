"""Module to write a basic detector GDML file

Attributes
----------
skeleton : str
    Multi-line string holding format for basic GDML file.
material_options : dict
    Dictionary of material names to their GDML definitions
    Materials must be named 'hunk_material' so that it can be
    referenced by the hunk volume
"""

import os

skeleton = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE gdml>
<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">
    <materials>
        <!-- We fill the non-hunk world with vaccuum -->
        <material name="vacuum" state="gas">
            <MEE unit="eV" value="19.2"/>
            <D unit="g/cm3" value="9.99999473841014e-09"/>
            <fraction n="1" ref="H"/>
        </material>
        <!-- We call whatever material we want to be in the hunk 'hunk_material' -->
        {hunk_material}
    </materials>

    <solids>
        <box name="world_box" x="{world_transverse}" y="{world_transverse}" z="{world_depth}" lunit="mm" />
        <box name="hunk_solid" x="{hunk_transverse}" y="{hunk_transverse}" z="{hunk_depth}" lunit="mm" />
    </solids>

    <structure>

      <!-- The EcalDarkBremFilter can be used if we name our volume in a specific way.
        That filter looks for the dark brem in all volumes that have 'volume' in the name
        AND one of 'W','Si','PCB', or 'CFMix'
        -->
        <volume name="hunk_W_volume">
            <materialref ref="hunk_material"/>
            <solidref ref="hunk_solid"/>
            <auxiliary auxtype="Region" auxvalue="CalorimeterRegion"/>
        </volume>

        <volume name="World">
            <materialref ref="vacuum"/>
            <solidref ref="world_box"/>
            
            <physvol name="hunk">
                <volumeref ref="hunk_W_volume"/>
                <position name="center" unit="mm" x="0" y="0" z="{hunk_depth}/2."/>
            </physvol> 

            <auxiliary auxtype="DetElem" auxvalue="Top"/>
        </volume>
    </structure>

    <userinfo>
        <auxiliary auxtype="DetectorVersion" auxvalue="999">
            <auxiliary auxtype="DetectorName" auxvalue="dark-brem-testing"/>
            <auxiliary auxtype="Description" 
                        auxvalue="A hunk of material to see how the dark brem process behaves."/>
        </auxiliary>

        <auxiliary auxtype="Region" auxvalue="CalorimeterRegion">
            <auxiliary auxtype="StoreTrajectories" auxvalue="false"/>
        </auxiliary>   
    </userinfo>  

    <setup name="Default" version="1.0">
        <world ref="World"/>
    </setup>
</gdml>
"""

material_options = {
    # LDMX and NA64 ECal absorber
    'tungsten' : """
        <isotope N="180" Z="74" name="W180">
            <atom unit="g/mole" value="179.947"/>
        </isotope>
        <isotope N="182" Z="74" name="W182">
            <atom unit="g/mole" value="181.948"/>
        </isotope>
        <isotope N="183" Z="74" name="W183">
            <atom unit="g/mole" value="182.95"/>
        </isotope>
        <isotope N="184" Z="74" name="W184">
            <atom unit="g/mole" value="183.951"/>
        </isotope>
        <isotope N="186" Z="74" name="W186">
            <atom unit="g/mole" value="185.954"/>
        </isotope>
        <element name="W">
            <fraction n="0.0012" ref="W180"/>
            <fraction n="0.265" ref="W182"/>
            <fraction n="0.1431" ref="W183"/>
            <fraction n="0.3064" ref="W184"/>
            <fraction n="0.2843" ref="W186"/>
        </element>
        <material name="hunk_material" state="solid">
            <T unit="K" value="293.15"/>
            <MEE unit="eV" value="727"/>
            <D unit="g/cm3" value="19.3"/>
            <fraction n="1" ref="W"/>
        </material>
        """,
    # lighter detection material
    'silicon' : """
        <isotope N="28" Z="14" name="Si280">
            <atom unit="g/mole" value="27.9769"/>
        </isotope>
        <isotope N="29" Z="14" name="Si290">
            <atom unit="g/mole" value="28.9765"/>
        </isotope>
        <isotope N="30" Z="14" name="Si300">
            <atom unit="g/mole" value="29.9738"/>
        </isotope>
        <element name="Si">
            <fraction n="0.922296077703922" ref="Si280"/>
            <fraction n="0.0468319531680468" ref="Si290"/>
            <fraction n="0.0308719691280309" ref="Si300"/>
        </element>
        <material name="hunk_material" state="solid">
            <T unit="K" value="293.15"/>
            <MEE unit="eV" value="173"/>
            <D unit="g/cm3" value="2.33"/>
            <fraction n="1" ref="Si"/>
        </material>
        """,
    # common test material?
    'iron' : """
      <isotope N="54" Z="26" name="Fe54">
         <atom unit="g/mole" value="53.9396"/>
      </isotope>
      <isotope N="56" Z="26" name="Fe56">
         <atom unit="g/mole" value="55.9349"/>
      </isotope>
      <isotope N="57" Z="26" name="Fe57">
         <atom unit="g/mole" value="56.9354"/>
      </isotope>
      <isotope N="58" Z="26" name="Fe58">
         <atom unit="g/mole" value="57.9333"/>
      </isotope>
      <element name="Fe">
         <fraction n="0.05845" ref="Fe54"/>
         <fraction n="0.91754" ref="Fe56"/>
         <fraction n="0.02119" ref="Fe57"/>
         <fraction n="0.00282" ref="Fe58"/>
      </element>
      <material name="hunk_material" state="solid">
        <T unit="K" value="293.15"/>
        <D unit="g/cm3" value="7.874"/>
        <fraction n="1" ref="Fe"/>
      </material>
      """,
    # BDX beam dump material
    'aluminum' : """
      <material Z="13" name="hunk_material" state="solid">
          <T unit="K" value="293.15"/>
          <MEE unit="eV" value="166"/>
          <D unit="g/cm3" value="2.699"/>
          <atom unit="g/mole" value="26.9815"/>
      </material>
      """,
    # CMS HCal absorber
    'brass' : """
      <element Z="29" name="Cu">
        <atom unit="g/mole" value="63.546"/>
        <D unit="g/cm3" value="8.96"/>
      </element>
      <element Z="30" name="Zn">
        <atom unit="g/mole" value="65.39"/>
        <D unit="g/cm3" value="7.112"/>
      </element>
      <material name="hunk_material" state="solid">
        <T unit="K" value="293.15"/>
        <D unit="g/cm3" value="8.53"/>
        <fraction n="0.7" ref="Cu"/>
        <fraction n="0.3" ref="Zn"/>
      </material>
      """,
    # LDMX and NA64 HCal absorber
    'steel' : """
      <isotope N="12" Z="6" name="C12">
         <atom unit="g/mole" value="12"/>
      </isotope>
      <isotope N="13" Z="6" name="C13">
         <atom unit="g/mole" value="13.0034"/>
      </isotope>
      <element name="C">
         <fraction n="0.9893" ref="C12"/>
         <fraction n="0.0107" ref="C13"/>
      </element>
      <isotope N="55" Z="25" name="Mn55">
         <atom unit="g/mole" value="54.938"/>
      </isotope>
      <element name="Mn">
         <fraction n="1" ref="Mn55"/>
      </element>
      <isotope N="54" Z="26" name="Fe54">
         <atom unit="g/mole" value="53.9396"/>
      </isotope>
      <isotope N="56" Z="26" name="Fe56">
         <atom unit="g/mole" value="55.9349"/>
      </isotope>
      <isotope N="57" Z="26" name="Fe57">
         <atom unit="g/mole" value="56.9354"/>
      </isotope>
      <isotope N="58" Z="26" name="Fe58">
         <atom unit="g/mole" value="57.9333"/>
      </isotope>
      <element name="Fe">
         <fraction n="0.05845" ref="Fe54"/>
         <fraction n="0.91754" ref="Fe56"/>
         <fraction n="0.02119" ref="Fe57"/>
         <fraction n="0.00282" ref="Fe58"/>
      </element>
      <material name="hunk_material" state="solid">
         <T unit="K" value="293.15"/>
         <D unit="g/cm3" value="7.87"/>
         <fraction n="0.9843" ref="Fe"/>
         <fraction n="0.014" ref="Mn"/>
         <fraction n="0.0017" ref="C"/>
      </material> 
      """,
      'lead' : """
      <material Z="82" name="hunk_material" state="solid">
          <T unit="K" value="293.15"/>
          <MEE unit="eV" value="823.0"/>
          <D unit="g/cm3" value="11.350"/>
          <atom unit="g/mole" value="207.2"/>
      </material>
      """
        }

def write(path, material, hunk_depth, hunk_transverse) :
    """Write the simple GDML detector file and return full path

    The world size is determined from the hunk size.
    The hunk is placed so that it's front face is at z = 0 
    so the primary particles can start just shy of z=0.

    Parameters
    ----------
    path : str
        path to new GDML file to write
    material : str
        name of material to use in hunk
        must be in the material_options dict
    hunk_depth : float
        Depth of hunk [mm]
    hunk_transverse : float
        Transverse width of hunk [mm]

    Returns
    -------
        str : full path to written GDML file

    Examples
    --------
    In your configuration script, if 'sim' is already defined
    to be a simulator, we can run the simulation within a hunk of
    tungsten using

        from LDMX.Detectors.write import write
        sim.detector = write('.write_detector.gdml','tungsten',500.,100.)
    """
    
    if material not in material_options :
        raise Exception(f'\'{material}\' is not one of the options: {list(material_options.keys())}')

    with open(path,'w') as gdml_f :
        gdml_f.write(skeleton.format(
            hunk_material = material_options[material],
            hunk_depth = hunk_depth,
            hunk_transverse = hunk_transverse,
            world_depth = 2*hunk_depth+1, 
            world_transverse = hunk_transverse+1
            ))
    return os.path.realpath(path)

if __name__ == '__main__' :
    import sys
    print(write('test_write.gdml',sys.argv[1], 500, 100))
