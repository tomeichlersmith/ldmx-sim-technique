# extra MG total xsec from set of MG/ME files

import pylhe
import os
import csv
import xml
el_d = '/local/cms/user/revering/madgraph/Repository/text/AN_doc/map0p1_el_W'
mu_d = '/local/cms/user/revering/madgraph/Repository/text/AN_doc/map_0p2_Cu'

def run_data(lhe_file) :
    """Get the beam energy [GeV] and cross section [pb]
    
    For some god-foresaken reason, the line in the MGGenerationInfo block
    labeled 'Integrated weight (pb)' and the entry in the init block that
    is the integrated weight by the schema definition are not the same always.
    
    We are choosing to just always use the MGGenerationInfo block one.
    """
    init_info = pylhe.readLHEInit(lhe_file)
    if len(init_info['procInfo']) == 0 :
        return None, None
    with open(lhe_file) as text :
        matches = [l for l in text.readlines() if 'Integrated' in l]
        if len(matches) != 1 :
            raise KeyError(f'More than one line in {lhe_file} with Integrated in it.')
    return init_info['initInfo']['energyA'], float(matches[0].split()[-1])

def extract_xsec(lhe_dir, output_file) :
    with open(output_file, 'w', newline='') as csvfile :
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['Energy [GeV]', 'Xsec [pb]'])
        for lhe_file in os.listdir(lhe_dir) :
            if not lhe_file.endswith('lhe') :
                continue
            try :
                beam, xsec = run_data(f'{lhe_dir}/{lhe_file}')
            except Exception :
                print(f'Skipping {lhe_file}')
                continue
            if beam is not None and xsec is not None :
                csvwriter.writerow([beam,xsec])

def main() :
    import sys
    extract_xsec(sys.argv[1], sys.argv[2])

if __name__ == '__main__' :
    main()
