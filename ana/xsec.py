"""Analyze the dark brem events passed to us on the command line."""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.ROOT)

from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

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
    inset_xlim = None, ymax = None,
    inset_width_height = (0.45,0.45),
    legend_overhead = False,
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
    
    if inset_xlim is not None :
        # we position the inset in the raw axes
        # and then mark it relative to the ratio axes we are zooming in on
        inset = plt.axes([0,0,1,1])
        # left edge, bottom edge, width, height
        inset.set_axes_locator(InsetPosition(raw, 
          [0.50,0.10,inset_width_height[0], inset_width_height[1]]))
        mark_inset(ratio, inset, loc1=2, loc2=4, fc='none', ec='0.5')

        selection = (mg_x > inset_xlim[0])&(mg_x < inset_xlim[1])
        inset.plot(mg_x[selection], [1 for x in mg_x[selection]], 
                   linewidth=0, label='MG/ME', 
                   #marker='o', markersize=10, markerfacecolor='none', 
                   markeredgecolor='tab:green')

    for name, data, style in others :
        y = data['Xsec [pb]']
        x = data['Energy [MeV]']/1000.
        raw.plot(x, y, label=name, linewidth=2, **style)
        data_interp = scipy.interpolate.interp1d(x, y)
        in_range = (mg_x >= x.min())&(mg_x <= x.max())#&(mg_x > 100)
        data_at_mge = [data_interp(e) for e in mg_x[in_range]]
        ratio_vals = data_at_mge/mg_y[in_range]
        ratio.plot(mg_x[in_range], ratio_vals,
            marker='.', markersize=15, linewidth=0, **style)

        if inset_xlim is not None :
            sl = in_range&(mg_x > inset_xlim[0])&(mg_x < inset_xlim[1])
            inset.plot(mg_x[sl], ratio_vals[sl], 
                marker='.', markersize=15, linewidth=0,
                label=name, **style)

    raw.set_ylabel('Total Cross Section / $\epsilon^2$ [pb]')
    raw.set_ylim(ymax=ymax)

    if legend_overhead :
        l = raw.legend(title=title, loc='lower center', bbox_to_anchor=(0.5,1))
        plt.setp(l.get_title(), multialignment='center')
    else :
        l = raw.legend(title=title, loc='upper left')
        plt.setp(l.get_title(), multialignment='left')

    if inset_xlim is not None :
        # give inset labels a white background so they cover the line
        #inset.set_xticklabels(inset.get_xticks(), backgroundcolor='w')
        pass

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

    xsec_plot('data/mg/el_xsec_with_lowE.csv', [
            ('G4DarkBreM Improved WW', 
              pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_iww.csv'),
              dict(color='tab:blue')),
          ],
         f'{arg.out_dir}/el_xsec.pdf',
        inset_xlim=(0,5), ymax=1.1e10,
        xlabel = 'Incident Electron Energy [GeV]',
        title = 'Electrons on Tungsten $m_{A\'} = 0.1$ GeV')

    xsec_plot('data/mg/mu_xsec_1GeV.csv', [
            ('G4DarkBreM Full WW', 
              pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_fullww.csv'),
              dict(color='tab:blue')),
          ],
        f'{arg.out_dir}/mu_xsec.pdf',
        inset_xlim=(0,100),
        xlabel = 'Incident Muon Energy [GeV]', ymax=5.5e6,
        title = 'Muons on Copper $m_{A\'} = 1.0$ GeV')

    dx, dy = plt.gcf().get_size_inches()
    plt.gcf().set_size_inches(dx, dy+4)

    xsec_plot('data/mg/el_xsec_with_lowE.csv', [
            ('G4DarkBreM Full WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_fullww.csv'),
              dict(color='tab:orange')),
            ('G4DarkBreM Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_iww.csv'),
              dict(color='tab:red')),
            ('G4DarkBreM Hyper-Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec_hiww.csv'),
              dict(color='tab:purple')),
          ],
         f'{arg.out_dir}/el_xsec_appendix.pdf',
        inset_xlim = (0,5), #legend_overhead=True,
        ymax = 13e9, inset_width_height=(0.45,0.25),
        xlabel = 'Incident Electron Energy [GeV]',
        title = 'Electrons on Tungsten $m_{A\'} = 0.1$ GeV')

    xsec_plot('data/mg/mu_xsec_1GeV.csv', [
            ('G4DarkBreM Full WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_fullww.csv'),
              dict(color='tab:orange')),
            ('G4DarkBreM Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_iww.csv'),
              dict(color='tab:red')),
            ('G4DarkBreM Hyper-Improved WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec_hiww.csv'),
              dict(color='tab:purple')),
          ],
        f'{arg.out_dir}/mu_xsec_appendix.pdf',
        #legend_overhead=True,
        inset_xlim=(0,100), ymax=9.5e6,
        inset_width_height=(0.45,0.25),
        xlabel = 'Incident Muon Energy [GeV]',
        title = 'Muons on Copper $m_{A\'} = 1.0$ GeV')

if __name__ == '__main__' :
    main()
