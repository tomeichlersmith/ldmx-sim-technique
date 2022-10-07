import uproot
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mplhep
import dark_brem_lhe

plt.style.use(mplhep.style.ROOT)

def read(f, tree_name = None) :
    df = None
    unit = 1.
    if f.endswith('.h5') :
        df = pd.read_hdf(f, key = tree_name)
    elif 'dblib' in f :
        if f.endswith('.lhe') :
            df = dark_brem_lhe.DarkBremEventFile(f).events
        else :
            df = dark_brem_lhe.DarkBremEventLibrary(f).events()
    else :
        if tree_name is None :
            raise KeyError('tree_name required for reading ROOT files')
        unit = 1000.
        with uproot.open(f'{f}:{tree_name}') as t :
            df = t.arrays(library='pd')

    df['recoil_px'] /= unit
    df['recoil_py'] /= unit
    df['recoil_pz'] /= unit
    df['recoil_energy'] /= unit
    df['recoil_pt'] = np.sqrt(df['recoil_px']**2 + df['recoil_py']**2)
    df['recoil_angle'] = np.arctan2(df['recoil_pt'],df['recoil_pz'])
    return df

def plt_efrac_angle(f, e_beam,
             angle_ax, kefrac_ax,
             angle_kw = {}, kefrac_kw = {}, 
             lepton_mass = None, tree_name = None) :
    """Pull out the energy and angle from the ROOT files
    Michael produced while studying the scaling technique

    Then do the calculations and binning, returning the binned
    data to save on memory space.
    """

    df = read(f, tree_name)

    weights = np.empty(len(df))
    weights.fill(1./len(df))

    ke = kefrac_ax.hist(df['recoil_energy']/e_beam, weights=weights, **kefrac_kw)
    ang = angle_ax.hist(df['recoil_angle'], weights=weights, **angle_kw)
    return ang, ke
    
def binned_hist(ax, np_hist, **kwargs) :
    vals, edges = np_hist
    ax.hist((edges[1:]+edges[:-1])/2, bins=edges, weights=vals, **kwargs)
    
def scaling_validation(file_packet,
                       ang_raw_ylim = None, ang_ratio_ylim = None, ke_raw_ylim = None, ke_ratio_ylim = None,
                       ke_legend_kw = dict(loc='lower right'), ang_legend_kw = dict(),
                       file_prefix = None) :
    (title, e_beam, (lepton, lepton_mass), mg, scaled) = file_packet

    (ke_raw, ke_ratio) = plt.figure('ke').subplots(ncols = 1, nrows = 2, 
                                                   sharex = 'col', gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ke_raw.set_ylabel('Fraction Below Outgoing Energy')
    if ke_raw_ylim is not None :
        ke_raw.set_ylim(ke_raw_ylim)
    ke_ratio.set_ylabel('Scaled / Unscaled')
    ke_ratio.set_xlabel(f'Outgoing {lepton} Energy Fraction')
    if ke_ratio_ylim is not None :
        ke_ratio.set_ylim(ke_ratio_ylim)

    (ang_raw, ang_ratio) = plt.figure('ang').subplots(ncols = 1, nrows = 2, 
                                                      sharex = 'col', gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ang_raw.set_yscale('log')
    ang_raw.set_ylabel('Normalized Rate')
    if ang_raw_ylim is not None :
        ang_raw.set_ylim(ang_raw_ylim)
    ang_ratio.set_ylabel('Scaled / Unscaled')
    ang_ratio.set_xlabel(f'Outgoing {lepton} Angle [rad]')
    if ang_ratio_ylim is not None :
        ang_ratio.set_ylim(ang_ratio_ylim)

    (ang_mg_vals, ang_mg_edges, ang_patches), (ke_mg_vals, ke_mg_edges, ke_patches) = plt_efrac_angle(mg,
        e_beam, ang_raw, ke_raw, lepton_mass = lepton_mass,
        angle_kw = dict(bins=25,range=(0.,3.14),label=f'Unscaled MG/ME at {e_beam} GeV',histtype='step',linewidth=2,color='black'), 
        kefrac_kw = dict(bins=50, range=(0.,1.), cumulative=True, label=f'Unscaled MG/ME at {e_beam} GeV',histtype='step',linewidth=2, color='black'))

    ke_ratio.axhline(1.,color='black')
    ang_ratio.axhline(1.,color='black')

    for name, f in scaled :
        (ang_vals, ang_edges, ang_patches), (ke_vals, ke_edges, ke_patches) = plt_efrac_angle(f, 
            e_beam, ang_raw, ke_raw, tree_name = 'dbint',
            angle_kw = dict(bins=25,range=(0.,3.14), label=name, histtype='step', linewidth=2), 
            kefrac_kw = dict(bins=50, range=(0.,1.), cumulative=True, label=name, histtype='step',linewidth=2))
            
        binned_hist(ke_ratio, (np.divide(ke_vals,ke_mg_vals), ke_edges), histtype='step', linewidth=2, label=name)
        binned_hist(ang_ratio, (np.divide(ang_vals,ang_mg_vals), ang_edges), histtype='step', linewidth=2, label=name)

    ke_raw.legend(title=title, **ke_legend_kw)
    ang_raw.legend(title=title, **ang_legend_kw)

    if file_prefix is not None :
        for fn in ['ke','ang'] :
            plt.figure(fn).savefig(f'{file_prefix}_{fn}.pdf')
            plt.figure(fn).clf()

def plt_efrac_pt(f, e_beam,
                 pt_ax, efrac_ax,
                 pt_kw = {}, efrac_kw = {}, 
                 lepton_mass = None, tree_name = None) :
    """Pull out the energy and angle from the ROOT files
    Michael produced while studying the scaling technique

    Then do the calculations and binning, returning the binned
    data to save on memory space.
    """

    df = read(f, tree_name)

    weights = np.empty(len(df))
    weights.fill(1./len(df))

    efrac = efrac_ax.hist(df['recoil_energy']/e_beam, weights = weights, **efrac_kw)
    pt = pt_ax.hist(df['recoil_pt'], weights=weights, **pt_kw)
    return pt, efrac

def efrac_pt_rates(data_packet,
                   efrac_legend_kw = dict(),
                   pt_legend_kw = dict(),
                   file_prefix = None) :
    (title, e_beam, data) = data_packet
    efrac_ax = plt.figure('efrac').subplots()
    pt_ax    = plt.figure('pt').subplots()
    
    efrac_ax.set_yscale('log')
    efrac_ax.set_ylabel('Normalized Rate')
    efrac_ax.set_xlabel('Outgoing Electron Energy Fraction')
    
    pt_ax.set_yscale('log')
    pt_ax.set_ylabel('Normalized Rate')
    pt_ax.set_xlabel('Outgoing Electron Transverse Momentum [GeV]')
    
    for name, f in data :
        plt_efrac_pt(f, e_beam,
                     pt_ax, efrac_ax,
                     pt_kw = dict(histtype='step',linewidth=2,range=(0,1.1),bins=50,label=name),
                     efrac_kw = dict(histtype='step',linewidth=2,range=(0,1),bins=50,label=name))
    
    efrac_ax.legend(title = title, **efrac_legend_kw)
    pt_ax.legend(title = title, **pt_legend_kw)
    
    if file_prefix is not None :
        for fn in ['efrac','pt'] :
            plt.figure(fn).savefig(f'{file_prefix}_{fn}.pdf')
            plt.figure(fn).clf()
    
def main() :
    electron = ('Electron', 0.000511)
    muon     = ('Muon'    , 0.105   )

    el_W = ('Electrons on Tungsten\n$m_{A\'} = 0.1$ GeV', 4., electron,
        'dblib/scaling/electron_tungsten_mA_0.1_E_4.0.h5',
        [
         ('Scaled from 4.2 GeV', 'data/scaling/el_W_4p2to4p0.root'),
         ('Scaled from 4.4 GeV', 'data/scaling/el_W_4p4to4p0.root'),
         ('Scaled from 4.8 GeV', 'data/scaling/el_W_4p8to4p0.root'),
         ('Scaled from 6.0 GeV', 'data/scaling/el_W_6p0to4p0.root'),
        ])

    mu_Pb = ('Muons on Lead\n$m_{A\'} = 1.0$ GeV', 100., muon, 
        'dblib/scaling/muon_lead_mA_1.0_E_100.h5',
       [('Scaled from 111 GeV', 'data/scaling/mu_Pb_111to100.root'),
        ('Scaled from 125 GeV', 'data/scaling/mu_Pb_125to100.root'),
        ('Scaled from 150 GeV', 'data/scaling/mu_Pb_150to100.root'),
        ('Scaled from 200 GeV', 'data/scaling/mu_Pb_200to100.root')
       ])
    
    scaling_validation(el_W,
                       ke_ratio_ylim=(0.975,1.025+0.003), 
                       ang_ratio_ylim=(0.8,1.25), ang_raw_ylim=(4e-4,4.),
                       ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
                       file_prefix='el_W')
    
    scaling_validation(mu_Pb,
                       ke_ratio_ylim=(0.8,1.25),
                       ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
                       file_prefix='mu_Pb')
    
    efrac_pt_rates(('4 GeV Electrons\n$m_{A\'} = 0.1$ GeV', 4.0,
        [
        ('W Target, Z = 74', 'dblib/scaling/electron_tungsten_mA_0.1_E_4.0.h5'),
        ('Cu Target, Z = 29', 'dblib/scaling/electron_copper_mA_0.1_E_4.0.h5'),
        ('Pb Target, Z = 82', 'dblib/scaling/electron_lead_mA_0.1_E_4.0.h5'),
        ]),
        file_prefix='el_material_comp')
    
    efrac_pt_rates(('100 GeV Muons\n$m_{A\'} = 1.0$ GeV', 100.0,
        [
        ('W Target, Z = 74', 'dblib/scaling/muon_tungsten_mA_1.0_E_100.h5'),
        ('Cu Target, Z = 29', 'dblib/scaling/muon_copper_mA_1.0_E_100.h5'),
        ('Pb Target, Z = 82', 'dblib/scaling/muon_lead_mA_1.0_E_100.h5'),
        ]),
        file_prefix='mu_material_comp')
    
if __name__ == '__main__' :
    main()
