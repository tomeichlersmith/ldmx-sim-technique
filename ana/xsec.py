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
import pickle as pkl

def parse_mg(csv_path) :
    df = pd.read_csv(csv_path) \
        .groupby('Energy [GeV]') \
        .apply(lambda samples : samples[samples >= samples.median() - 2*samples.std()].mean()) \
        .drop(columns='Energy [GeV]') \
        .reset_index()
    df['Xsec [pb]'] *= (127.9/137)**3
    df['Energy [MeV]'] = df['Energy [GeV]']*1000
    return df

def xsec_plot(mg_csv, others, file_name, 
    xlabel = 'Incident Lepton Energy [GeV]', title = None) :
    (raw, ratio) = plt.gcf().subplots(ncols = 1, nrows = 2, 
        sharex = 'col', gridspec_kw=dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace=0)

    mg = parse_mg(mg_csv)

    mg_x = mg['Energy [GeV]']
    mg_y = mg['Xsec [pb]']

    raw.plot(mg_x, mg_y, linewidth=0, label='MG/ME',
        marker='o', markersize=10, markerfacecolor='none', markeredgecolor='tab:green')
    ratio.plot(mg_x, [1. for x in mg_x], marker='.', markersize=0, linewidth=0, color='tab:green')

    for name, data, style in others :
        y = data['Xsec [pb]']
        x = data['Energy [MeV]']/1000.
        raw.plot(x, y, label=name, linewidth=2, **style)
        data_interp = scipy.interpolate.interp1d(x, y)
        in_range = (mg_x >= x.min())&(mg_x <= x.max())#&(mg_x > 100)
        data_at_mge = [data_interp(e) for e in mg_x[in_range]]
        ratio.plot(mg_x[in_range], data_at_mge/mg_y[in_range], 
            marker='.', markersize=15, linewidth=0, **style)

    #ratio.axhline(0.9, color='gray')

    raw.set_ylabel('Total Cross Section / $\epsilon^2$ [pb]')
    l = raw.legend(title=title)
    plt.setp(l.get_title(), multialignment='right')

    ratio.set_ylabel('Ratio to MG/ME')
    ratio.set_xlabel(xlabel)
    plt.savefig(file_name, bbox_inches='tight')
    with open(file_name.replace('pdf','pkl'),'wb') as f :
        pkl.dump(plt.gcf(), f)
    plt.clf()

def main() :
    import argparse
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir',help='Directory data is in', default='data/dev')
    parser.add_argument('--out_dir',help='Directory to put plots (Default: data_dir)')
    
    arg = parser.parse_args()
    
    # make sure output directory exists
    if arg.out_dir is None :
        arg.out_dir = arg.data_dir
    os.makedirs(arg.out_dir, exist_ok=True)

    xsec_plot('data/mg/el_xsec.csv', [
            ('G4DarkBreM Improved WW', 
              pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_iww.csv'),
              dict(color='tab:blue')),
          ],
         f'{arg.out_dir}/el_xsec.pdf',
        xlabel = 'Incident Electron Energy [GeV]',
        title = '$m_{A\'} = 0.1$ GeV\nElectrons on Tungsten')

    xsec_plot('data/mg/mu_xsec_1GeV.csv', [
            ('G4DarkBreM Full WW', 
              pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_fullww.csv'),
              dict(color='tab:blue')),
          ],
        f'{arg.out_dir}/mu_xsec.pdf',
        xlabel = 'Incident Muon Energy [GeV]',
        title = '$m_{A\'} = 1.0$ GeV\nMuons on Copper')

    xsec_plot('data/mg/el_xsec.csv', [
            ('G4DarkBreM Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_iww.csv'),
              dict(color='tab:blue')),
            ('G4DarkBreM Hyper-Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_hiww.csv'),
              dict(color='tab:orange')),
            ('G4DarkBreM Full WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_fullww.csv'),
              dict(color='tab:red')),
          ],
         f'{arg.out_dir}/el_xsec_appendix.pdf',
        xlabel = 'Incident Electron Energy [GeV]',
        title = '$m_{A\'} = 0.1$ GeV\nElectrons on Tungsten')

    xsec_plot('data/mg/mu_xsec_1GeV.csv', [
            ('G4DarkBreM Full WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_fullww.csv'),
              dict(color='tab:blue')),
            ('G4DarkBreM Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_iww.csv'),
              dict(color='tab:orange')),
            ('G4DarkBreM Hyper-Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_hiww.csv'),
              dict(color='tab:red')),
          ],
        f'{arg.out_dir}/mu_xsec_appendix.pdf',
        xlabel = 'Incident Muon Energy [GeV]',
        title = '$m_{A\'} = 1.0$ GeV\nMuons on Copper')

if __name__ == '__main__' :
    main()
