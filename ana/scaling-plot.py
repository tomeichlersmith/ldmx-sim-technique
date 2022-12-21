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
    elif f.endswith('.csv') :
        unit = 1000.
        df = pd.read_csv(f)
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
    df.drop(df[df.recoil_energy <= 0.0].index, inplace=True)
    return df
    
def plt_bins(ax, np_hist, **kwargs) :
    vals, edges = np_hist
    return ax.hist((edges[1:]+edges[:-1])/2, bins=edges, weights=vals, **kwargs)
    
def binit(data, **kwargs) :
    cumulative = kwargs.get('cumulative',False)
    if cumulative :
        del kwargs['cumulative']
    v, b = np.histogram(data, **kwargs)
    if cumulative :
        v = np.cumsum(v)
    # total number of entries
    total = v[-1]
    # error in cumulative bins is error
    #   in count above that bin
    e = np.sqrt(total-v)/v
    # normalize
    v = v / total
    print(v)
    print(e)
    return v, b, e

def efrac_cumulative_ratios(material_packs) :
    binned = {}
    e_beam = None
    for material, e, mg, scaled in material_packs :
        if e_beam is None :
            e_beam = e
        if e_beam != e :
            raise KeyError('Material pack with different beam energies.')
            
        mgdf = read(mg)
        mg_vals, mg_bins, mg_errs = binit(mgdf['recoil_energy']/e_beam,
                                          range=(0,1), bins=50, cumulative=True)
        ratios = []
        for name, f in scaled :
            df = read(f, 'dbint')
            vals, bins, errs = binit(df['recoil_energy']/e_beam,
                                     range=(0,1), bins=50, cumulative=True)
            ratio = vals/mg_vals
            #error = ratio*np.sqrt(errs**2 + mg_errs**2)
            #error = np.sqrt(errs**2 + mg_errs**2)
            error = mg_errs
            ratios.append((name, ratio, error))
        binned[material] = (mg_bins, ratios)
    return e_beam, binned

def plt_cumulative_ratios(lepton, e_beam, binned, file_prefix, 
                          ax_kwargs = None, ap = 0.1,
                          material_label_y = 0.25) :
    fig = plt.figure(figsize=(8,11))
    axes = fig.subplots(ncols = 1, nrows = len(binned), sharex = 'col')
    if ax_kwargs is None :
        ax_kwargs = [dict() for i in range(len(binned))]
    for i, (ax, ax_kw, hists) in enumerate(zip(axes, ax_kwargs, binned.items())) :
        ax.set_ylim(*ax_kw.get('ylim',(0.965,1.025)))
        material, (bins, vals) = hists
        ax.text(0.9, material_label_y, material, ha='right', va='top', transform=ax.transAxes)
        ax.axhline(1,color='black')
        for name, v, er in vals :
            lines = ax.plot((bins[:-1]+bins[1:])/2, v, 
                            linewidth=2, drawstyle='steps-mid',
                            label=name if i == 0 else '_nolegend_')
            len(lines)
            print(er)
            ax.bar(bins[:-1], er, align='edge', 
                   width = (bins[1:]-bins[:-1]), 
                   bottom = 1, color='lightgray')
            ax.bar(bins[:-1], -1*er, align='edge', 
                   width = (bins[1:]-bins[:-1]), 
                   bottom = 1, color='lightgray')
        
        if i == 0 :
            ax.set_ylabel('Ratio Scaled / Unscaled Predicted')
        
        if i == len(binned)-1 :
            ax.set_xlabel(f'Maximum Outgoing {lepton[0]} Energy Fraction')
            
    fig.legend(
        title='$m_{A\'} = '+str(ap)+'$ GeV\n'+f'Scaled to\n{e_beam} GeV from',
        bbox_to_anchor=(0.9,0.5), loc='center left')
    fig.savefig(f'{file_prefix}_efrac_cumulative_ratio.pdf', bbox_inches='tight')
    fig.clf()

def plt_efrac_angle(f, e_beam,
             angle_ax, kefrac_cumu_ax, kefrac_diff_ax,
             angle_kw = {}, kefrac_kw = {}, 
             lepton_mass = None, tree_name = None) :
    """Pull out the energy and angle from the ROOT files
    Michael produced while studying the scaling technique

    Then do the calculations and binning, returning the binned
    data to save on memory space.
    """

    df = read(f, tree_name)

    weights = np.ones(len(df))/len(df)

    ke = kefrac_diff_ax.hist(df['recoil_energy']/e_beam, weights = weights, **kefrac_kw)
    ke_cumu = kefrac_cumu_ax.hist(df['recoil_energy']/e_beam, weights = weights, **kefrac_kw, cumulative=True)
    ang = angle_ax.hist(df['recoil_angle'], weights = weights, **angle_kw)
    return ang, ke, ke_cumu

def scaling_validation(file_packet, apmass, lepton_packet,
                       ang_raw_ylim = None, ang_ratio_ylim = None, 
                       ke_ratio_ylim = None,
                       ke_legend_kw = dict(loc='lower right'), ang_legend_kw = dict(),
                       file_prefix = None) :
  
    kefrac_bins = np.arange(0,1,1/50) #np.concatenate((np.arange(0,10,1),np.arange(10,100,10),np.arange(100,4000,100)))/4000.
    (material, e_beam, mg, scaled) = file_packet
    (lepton, lepton_mass) = lepton_packet

    (ke_cumu_raw, ke_cumu_ratio) = plt.figure('ke_cumu').subplots(ncols = 1, nrows = 2, 
                                                   sharex = 'col', 
                                                   gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ke_cumu_raw.set_ylabel('Fraction Below Outgoing Energy')
    ke_cumu_ratio.set_ylabel('Scaled / Unscaled')
    ke_cumu_ratio.set_xlabel(f'Outgoing {lepton} Energy Fraction')
    #ke_cumu_raw.set_xscale('log')
    if ke_ratio_ylim is not None :
        ke_cumu_ratio.set_ylim(ke_ratio_ylim)


    (ke_diff_raw, ke_diff_ratio) = plt.figure('ke_diff').subplots(ncols = 1, nrows = 2, 
                                                   sharex = 'col', 
                                                   gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ke_diff_raw.set_ylabel('Normalized Rate')
    ke_diff_ratio.set_ylabel('Scaled / Unscaled')
    ke_diff_ratio.set_xlabel(f'Outgoing {lepton} Energy Fraction')
    #ke_diff_raw.set_xscale('log')
    if ke_ratio_ylim is not None :
        ke_diff_ratio.set_ylim(ke_ratio_ylim)

    (ang_raw, ang_ratio) = plt.figure('ang').subplots(ncols = 1, nrows = 2, 
                                                      sharex = 'col', 
                                                      gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ang_raw.set_yscale('log')
    ang_raw.set_ylabel('Normalized Rate')
    if ang_raw_ylim is not None :
        ang_raw.set_ylim(ang_raw_ylim)
    ang_ratio.set_ylabel('Scaled / Unscaled')
    ang_ratio.set_xlabel(f'Outgoing {lepton} Angle [rad]')
    if ang_ratio_ylim is not None :
        ang_ratio.set_ylim(ang_ratio_ylim)

    (ang_mg_vals, ang_mg_edges, ang_patches), (ke_mg_vals, ke_mg_edges, ke_patches), (ke_mg_cumu_vals, ke_mg_cumu_edges, ke_mg_cumu_patches) = plt_efrac_angle(mg,
        e_beam, ang_raw, ke_cumu_raw, ke_diff_raw,
        lepton_mass = lepton_mass,
        angle_kw = dict(bins=25,range=(0.,3.14),
          label=f'Unscaled MG/ME at {e_beam} GeV',
          histtype='step',linewidth=2,color='black'), 
        kefrac_kw = dict(bins=kefrac_bins, density=True,
          label=f'Unscaled MG/ME at {e_beam} GeV',
          histtype='step',linewidth=2, color='black'))

    ke_cumu_ratio.axhline(1.,color='black')
    ke_diff_ratio.axhline(1.,color='black')
    ang_ratio.axhline(1.,color='black')

    for name, f in scaled :
        (ang_vals, ang_edges, ang_patches), (ke_vals, ke_edges, ke_patches), (ke_cumu_vals, ke_cumu_edges, ke_cumu_patches) = plt_efrac_angle(f,
            e_beam, ang_raw, ke_cumu_raw, ke_diff_raw,
            tree_name = 'dbint',
            angle_kw = dict(bins=25,range=(0.,3.14), 
              label=f'Scaled from {name}', histtype='step', linewidth=2), 
            kefrac_kw = dict(bins=kefrac_bins, density=True,
              label=f'Scaled from {name}', histtype='step',linewidth=2))
            
        plt_bins(ke_diff_ratio, (np.divide(ke_vals,ke_mg_vals), ke_edges), 
            histtype='step', linewidth=2, label=name)
        plt_bins(ke_cumu_ratio, (np.divide(ke_cumu_vals,ke_mg_cumu_vals), ke_cumu_edges), 
            histtype='step', linewidth=2, label=name)
        plt_bins(ang_ratio, (np.divide(ang_vals,ang_mg_vals), ang_edges), 
            histtype='step', linewidth=2, label=name)

    title = f'{lepton}s on {material}'+'\n$m_{A\'} = '+str(apmass)+'$ GeV'
    ke_cumu_raw.legend(title=title, **ke_legend_kw)
    ke_diff_raw.legend(title=title, **ke_legend_kw)
    ang_raw.legend(title=title, **ang_legend_kw)

    if file_prefix is not None :
        for fn in ['ke_cumu','ke_diff','ang'] :
            plt.figure(fn).savefig(f'{file_prefix}_{fn}.pdf', bbox_inches='tight')
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
    (title, data) = data_packet
    efrac_ax = plt.figure('efrac', figsize=(11,8)).subplots()
    pt_ax    = plt.figure('pt', figsize=(11,8)).subplots()
    
    efrac_ax.set_yscale('log')
    efrac_ax.set_ylabel('Normalized Rate')
    efrac_ax.set_xlabel('Outgoing Electron Energy Fraction')
    
    pt_ax.set_yscale('log')
    pt_ax.set_ylabel('Normalized Rate')
    pt_ax.set_xlabel('Outgoing Electron Transverse Momentum [GeV]')
    
    for name, e_beam, f in data :
        plt_efrac_pt(f, e_beam,
                     pt_ax, efrac_ax,
                     pt_kw = dict(histtype='step',linewidth=2,range=(0,1.1),bins=50,label=name),
                     efrac_kw = dict(histtype='step',linewidth=2,range=(0,1),bins=50,label=name))
    
    efrac_ax.legend(title = title, **efrac_legend_kw)
    pt_ax.legend(title = title, **pt_legend_kw)
    
    if file_prefix is not None :
        for fn in ['efrac','pt'] :
            plt.figure(fn).savefig(f'{file_prefix}_{fn}.pdf',bbox_inches='tight')
            plt.figure(fn).clf()
    
def main() :
    electron = ('Electron', 0.000511)
    muon     = ('Muon'    , 0.105   )

    import os
    os.makedirs('data/scaling/plots/', exist_ok=True)

    el_W = ('Tungsten', 4.,
        'dblib/scaling/electron_tungsten_mA_0.1_E_4.0.h5',
        [
         ('4.2 GeV', 'data/scaling/electron_tungsten_4.2_to_4.0.csv'),
         ('4.4 GeV', 'data/scaling/electron_tungsten_4.4_to_4.0.csv'),
         ('4.8 GeV', 'data/scaling/electron_tungsten_4.8_to_4.0.csv'),
         ('6.0 GeV', 'data/scaling/electron_tungsten_6.0_to_4.0.csv'),
        ])
    scaling_validation(el_W, 0.1, electron,
        ke_ratio_ylim=(0.975,1.025+0.003), 
        ang_ratio_ylim=(0.8,1.25), ang_raw_ylim=(4e-4,4.),
        ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
        file_prefix='data/scaling/plots/electron_tungsten')

    el_Pb = ('Lead', 4.,
        'dblib/scaling/electron_lead_mA_0.1_E_4.0.h5',
        [
         ('4.2 GeV', 'data/scaling/electron_lead_4.2_to_4.0.csv'),
         ('4.4 GeV', 'data/scaling/electron_lead_4.4_to_4.0.csv'),
         ('4.8 GeV', 'data/scaling/electron_lead_4.8_to_4.0.csv'),
         ('6.0 GeV', 'data/scaling/electron_lead_6.0_to_4.0.csv'),
        ])
    scaling_validation(el_Pb, 0.1, electron,
        ke_ratio_ylim=(0.965,1.015+0.003),
        ang_ratio_ylim=(0.8,1.25), ang_raw_ylim=(4e-4,4.),
        ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
        file_prefix='data/scaling/plots/electron_lead')

    el_Cu = ('Copper', 4.,
        'dblib/scaling/electron_copper_mA_0.1_E_4.0.h5',
        [
         ('4.2 GeV', 'data/scaling/electron_copper_4.2_to_4.0.csv'),
         ('4.4 GeV', 'data/scaling/electron_copper_4.4_to_4.0.csv'),
         ('4.8 GeV', 'data/scaling/electron_copper_4.8_to_4.0.csv'),
         ('6.0 GeV', 'data/scaling/electron_copper_6.0_to_4.0.csv'),
        ])
    scaling_validation(el_Cu, 0.1, electron,
        ke_ratio_ylim=(0.955,1.015+0.003),
        ang_ratio_ylim=(0.8,1.25), ang_raw_ylim=(4e-4,4.),
        ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
        file_prefix='data/scaling/plots/electron_copper')

    e_beam, binned = efrac_cumulative_ratios([el_Cu, el_W, el_Pb])
    plt_cumulative_ratios(electron, e_beam, binned, 'electron',
                          ax_kwargs=[
                            {'ylim':(0.955,1.01)},
                            {'ylim':(0.965,1.01)},
                            {'ylim':(0.965,1.01)}
                            ])

    efrac_pt_rates(('4 GeV Electrons\n$m_{A\'} = 0.1$ GeV',
        [
        ('Cu Target, Z = 29', 4.0, 'dblib/scaling/electron_copper_mA_0.1_E_4.0.h5'),
        ('W Target, Z = 74', 4.0, 'dblib/scaling/electron_tungsten_mA_0.1_E_4.0.h5'),
        ('Pb Target, Z = 82', 4.0, 'dblib/scaling/electron_lead_mA_0.1_E_4.0.h5'),
        ]),
        file_prefix='data/scaling/plots/electron_material_comp',
        efrac_legend_kw = dict(bbox_to_anchor=(0.05,0),loc='lower left'))

    efrac_pt_rates(('Electrons on Tungsten\n$m_{A\'} = 0.1$ GeV',
      [
        ('4.0 GeV', 4.0, 'dblib/scaling/electron_tungsten_mA_0.1_E_4.0.h5'),
        ('4.4 GeV', 4.4, 'dblib/scaling/electron_tungsten_mA_0.1_E_4.4.h5'),
        ('6.0 GeV', 6.0, 'dblib/scaling/electron_tungsten_mA_0.1_E_6.0.h5'),
        ('8.0 GeV', 8.0, 'dblib/scaling/electron_tungsten_mA_0.1_E_8.0.h5'),
      ]),
      efrac_legend_kw = dict(loc='lower left', bbox_to_anchor=(0.05,0.)),
      file_prefix='data/scaling/plots/electron_energy_comp')

    mu_W = ('Tungsten', 100.,
        'dblib/scaling/muon_tungsten_mA_1.0_E_100.h5',
       [('111 GeV', 'data/scaling/muon_tungsten_111_to_100.root'),
        ('125 GeV', 'data/scaling/muon_tungsten_125_to_100.root'),
        ('150 GeV', 'data/scaling/muon_tungsten_150_to_100.root'),
        ('200 GeV', 'data/scaling/muon_tungsten_200_to_100.root')
       ])
    scaling_validation(mu_W, 1.0, muon,
       ke_ratio_ylim=(0.95,1.25+0.003),
       ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
       file_prefix='data/scaling/plots/muon_tungsten')

    mu_Pb = ('Lead', 100.,
        'dblib/scaling/muon_lead_mA_1.0_E_100.h5',
       [('111 GeV', 'data/scaling/muon_lead_111_to_100.root'),
        ('125 GeV', 'data/scaling/muon_lead_125_to_100.root'),
        ('150 GeV', 'data/scaling/muon_lead_150_to_100.root'),
        ('200 GeV', 'data/scaling/muon_lead_200_to_100.root')
       ])
    scaling_validation(mu_Pb, 1.0, muon,
       ke_ratio_ylim=(0.95,1.25+0.003),
       ke_legend_kw=dict(bbox_to_anchor=(0.95,0), loc='lower right'),
       file_prefix='data/scaling/plots/muon_lead')

    mu_Cu = ('Copper', 100.,
        'dblib/scaling/muon_copper_mA_1.0_E_100.h5',
       [('111 GeV', 'data/scaling/muon_copper_111_to_100.root'),
        ('125 GeV', 'data/scaling/muon_copper_125_to_100.root'),
        ('150 GeV', 'data/scaling/muon_copper_150_to_100.root'),
        ('200 GeV', 'data/scaling/muon_copper_200_to_100.root')
       ])
    scaling_validation(mu_Cu, 1.0, muon,
       ke_ratio_ylim=(0.95,1.25+0.003),
       ke_legend_kw=dict(bbox_to_anchor=(0.95,0), loc='lower right'),
       file_prefix='data/scaling/plots/muon_copper')

    e_beam, binned = efrac_cumulative_ratios([mu_Cu, mu_W, mu_Pb])
    plt_cumulative_ratios(muon, e_beam, binned, 'muon', ap = 1.0,
                          material_label_y = 0.75,
                          ax_kwargs=[
                            {'ylim':(0.950,1.25)},
                            {'ylim':(0.950,1.25)},
                            {'ylim':(0.950,1.25)}
                            ])
    
    efrac_pt_rates(('100 GeV Muons\n$m_{A\'} = 1.0$ GeV',
        [
        ('Cu Target, Z = 29', 100, 'dblib/scaling/muon_copper_mA_1.0_E_100.h5'),
        ('W Target, Z = 74', 100, 'dblib/scaling/muon_tungsten_mA_1.0_E_100.h5'),
        ('Pb Target, Z = 82', 100, 'dblib/scaling/muon_lead_mA_1.0_E_100.h5'),
        ]),
        file_prefix='data/scaling/plots/muon_material_comp',
        efrac_legend_kw=dict(bbox_to_anchor=(0.05,0),loc='lower left'))

    efrac_pt_rates(('Muons on Copper\n$m_{A\'} = 1.0$ GeV',
      [
        ('100 GeV', 100, 'dblib/scaling/muon_copper_mA_1.0_E_100.h5'),
        ('111 GeV', 111, 'dblib/scaling/muon_copper_mA_1.0_E_111.h5'),
        ('125 GeV', 125, 'dblib/scaling/muon_copper_mA_1.0_E_125.h5'),
        ('150 GeV', 150, 'dblib/scaling/muon_copper_mA_1.0_E_150.h5'),
        ('200 GeV', 200, 'dblib/scaling/muon_copper_mA_1.0_E_200.h5'),
      ]),
      file_prefix='data/scaling/plots/muon_energy_comp')
    
if __name__ == '__main__' :
    main()
