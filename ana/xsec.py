"""Analyze the dark brem events passed to us on the command line."""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.ROOT)

import uproot
import numpy as np
import pandas as pd
import dark_brem_lhe
import scipy

def xsec_plot(mg_csv, others, file_name, 
    xlabel = 'Incident Lepton Energy [GeV]', title = None) :
    (raw, ratio) = plt.gcf().subplots(ncols = 1, nrows = 2, 
        sharex = 'col', gridspec_kw=dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace=0)

    mg = pd.read_csv(mg_csv)
    if 'el' in mg_csv :
        # clean out "bad" MG runs from estimate
        mg = mg.groupby('Energy [GeV]').apply(
            lambda samples : samples[samples >= samples.median() - 2*samples.std()].mean()).drop(columns='Energy [GeV]').reset_index()

    mg_x = mg['Energy [GeV]']
    mg_y = mg['Xsec [pb]']*(127.9/137)**3

    for name, data in others :
        y = data['Xsec [pb]']
        x = data['Energy [MeV]']/1000.
        raw.plot(x, y, label=name, linewidth=2)
        data_interp = scipy.interpolate.interp1d(x, y)
        in_range = (mg_x >= x.min())&(mg_x <= x.max())
        data_at_mge = [data_interp(e) for e in mg_x[in_range]]
        ratio.plot(mg_x[in_range], data_at_mge/mg_y[in_range], marker='.', markersize=15, linewidth=0)

    raw.plot(mg_x, mg_y,marker='.', markersize=15, linewidth=0, label='MG/ME')
    ratio.plot(mg_x, [1. for x in mg_x], marker='.', markersize=0, linewidth=0)

    raw.set_ylabel('Total Cross Section / $\epsilon^2$ [pb]')
    l = raw.legend(title=title)
    plt.setp(l.get_title(), multialignment='right')

    ratio.set_ylabel('Ratio to MG/ME')
    ratio.set_xlabel(xlabel)
    plt.savefig(file_name, bbox_inches='tight')
    plt.clf()

def main() :
    import argparse
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir',help='Directory data is in')
    parser.add_argument('--out_dir',help='Directory to put plots (Default: data_dir)')
    
    arg = parser.parse_args()
    
    # make sure output directory exists
    if arg.out_dir is None :
        arg.out_dir = arg.data_dir
    os.makedirs(arg.out_dir, exist_ok=True)

    xsec_plot('data/mg/mu_xsec_0.2GeV.csv', [
            ('G4DarkBreM', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_0.2GeV.csv')),
            ('DMG4', pd.read_csv(f'data/dev/dmg4_mu_xsec_0.2GeV.csv')),
          ],
        f'{arg.out_dir}/mu_xsec_0.2GeV.pdf',
        xlabel = 'Incident Muon Energy [GeV]',
        title = '$m_{A\'} = 0.2$ GeV\nMuons on Copper')

    xsec_plot('data/mg/mu_xsec_1GeV.csv', [
            ('G4DarkBreM', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec.csv')),
            ('DMG4', pd.read_csv(f'{arg.data_dir}/dmg4_mu_xsec_1.0GeV.csv')),
          ],
        f'{arg.out_dir}/mu_xsec_1GeV.pdf',
        xlabel = 'Incident Muon Energy [GeV]',
        title = '$m_{A\'} = 1.0$ GeV\nMuons on Copper')


    xsec_plot('data/mg/el_xsec.csv', [
            ('G4DarkBreM', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec.csv')),
          ],
         f'{arg.out_dir}/el_xsec.pdf',
        xlabel = 'Incident Electron Energy [GeV]',
        title = '$m_{A\'} = 0.1$ GeV\nElectrons on Tungsten')

if __name__ == '__main__' :
    main()
