import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplhep

plt.style.use(mplhep.style.ROOT)

repository = '/local/cms/user/revering/dphoton/ldmx-sim-technique/michaelCode/plotting/repository/'

def plt_efrac_angle(root_file, tree_name, e_beam,
             angle_ax, kefrac_ax,
             angle_kw = {}, kefrac_kw = {}, lepton_mass = None) :
    """Pull out the energy and angle from the ROOT files
    Michael produced while studying the scaling technique

    Then do the calculations and binning, returning the binned
    data to save on memory space.
    """

    with uproot.open(f'{repository}{root_file}:{tree_name}') as t :
        data_array = t.arrays('IncidentParticle', entry_stop=1000000, library='pd')
        x = data_array['IncidentParticle']['fP']['fX']
        y = data_array['IncidentParticle']['fP']['fY']
        z = data_array['IncidentParticle']['fP']['fZ']
        pt = np.sqrt(x*x + y*y)
        
        lepton_E = data_array['IncidentParticle']['fE']
        if lepton_mass is None :
            kefrac = lepton_E/e_beam
        else :
            kefrac = (lepton_E - lepton_mass)/e_beam
        angle  = np.arctan2(pt,z)
        
        weights = np.empty(len(kefrac))
        weights.fill(1./len(kefrac))

        ke = kefrac_ax.hist(kefrac, weights=weights, **kefrac_kw)
        ang = angle_ax.hist(angle, weights=weights, **angle_kw)
        return ang, ke
    
def binned_hist(ax, np_hist, **kwargs) :
    vals, edges = np_hist
    ax.hist((edges[1:]+edges[:-1])/2, bins=edges, weights=vals, **kwargs)
    
def scaling_validation(file_packet, lepton, title, e_beam, 
                       ang_raw_ylim = None, ang_ratio_ylim = None, ke_raw_ylim = None, ke_ratio_ylim = None,
                       ke_legend_kw = dict(loc='lower right'), ang_legend_kw = dict(),
                       file_prefix = None) :
    (mg, scaled) = file_packet

    (ke_raw, ke_ratio) = plt.figure('ke').subplots(ncols = 1, nrows = 2, 
                                                   sharex = 'col', 
                                                   gridspec_kw = dict(
                                                       height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ke_raw.set_ylabel('Fraction Below Outgoing Energy')
    if ke_raw_ylim is not None :
        ke_raw.set_ylim(ke_raw_ylim)
    ke_ratio.set_ylabel('Scaled / Unscaled')
    ke_ratio.set_xlabel(f'Outgoing {lepton} Energy Fraction')
    if ke_ratio_ylim is not None :
        ke_ratio.set_ylim(ke_ratio_ylim)

    (ang_raw, ang_ratio) = plt.figure('ang').subplots(ncols = 1, nrows = 2, sharex = 'col', gridspec_kw = dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace = 0)
    ang_raw.set_yscale('log')
    ang_raw.set_ylabel('Normalized Rate')
    if ang_raw_ylim is not None :
        ang_raw.set_ylim(ang_raw_ylim)
    ang_ratio.set_ylabel('Scaled / Unscaled')
    ang_ratio.set_xlabel(f'Outgoing {lepton} Angle [rad]')
    if ang_ratio_ylim is not None :
        ang_ratio.set_ylim(ang_ratio_ylim)

    (ang_mg_vals, ang_mg_edges, ang_patches), (ke_mg_vals, ke_mg_edges, ke_patches) = plt_efrac_angle(mg, 'Events', e_beam, ang_raw, ke_raw,
                             angle_kw = dict(bins=25,range=(0.,3.14),label=f'Unscaled MG/ME at {e_beam} GeV',histtype='step',linewidth=2,color='black'), 
                             kefrac_kw = dict(bins=50, range=(0.,1.), cumulative=True, label=f'Unscaled MG/ME at {e_beam} GeV',histtype='step',linewidth=2, color='black'))

    ke_ratio.axhline(1.,color='black')
    ang_ratio.axhline(1.,color='black')

    for name, f in scaled :
        (ang_vals, ang_edges, ang_patches), (ke_vals, ke_edges, ke_patches) = plt_efrac_angle(f, 'forward_only', e_beam, ang_raw, ke_raw,
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

def plt_efrac_pt(root_file, tree_name, e_beam,
                 pt_ax, efrac_ax,
                 pt_kw = {}, efrac_kw = {}, lepton_mass = None) :
    """Pull out the energy and angle from the ROOT files
    Michael produced while studying the scaling technique

    Then do the calculations and binning, returning the binned
    data to save on memory space.
    """

    with uproot.open(f'{repository}{root_file}:{tree_name}') as t :
        data_array = t.arrays('IncidentParticle', entry_stop=1000000, library='pd')
        x = data_array['IncidentParticle']['fP']['fX']
        y = data_array['IncidentParticle']['fP']['fY']
        z = data_array['IncidentParticle']['fP']['fZ']
        pt = np.sqrt(x*x + y*y)
        
        lepton_E = data_array['IncidentParticle']['fE']
        if lepton_mass is None :
            efrac = lepton_E/e_beam
        else :
            efrac = (lepton_E - lepton_mass)/e_beam
        
        weights = np.empty(len(efrac))
        weights.fill(1./len(efrac))

        efrac = efrac_ax.hist(efrac, weights=weights, **efrac_kw)
        pt = pt_ax.hist(pt, weights=weights, **pt_kw)
        return pt, efrac

def efrac_pt_rates(data, title, e_beam, 
                   efrac_legend_kw = dict(),
                   pt_legend_kw = dict(),
                   file_prefix = None) :
    efrac_ax = plt.figure('efrac').subplots()
    pt_ax    = plt.figure('pt').subplots()
    
    efrac_ax.set_yscale('log')
    efrac_ax.set_ylabel('Normalized Rate')
    efrac_ax.set_xlabel('Outgoing Electron Energy Fraction')
    
    pt_ax.set_yscale('log')
    pt_ax.set_ylabel('Normalized Rate')
    pt_ax.set_xlabel('Outgoing Electron Transverse Momentum [GeV]')
    
    for name, file in data :
        plt_efrac_pt(file, 'Events', e_beam,
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
    import sys
    if len(sys.argv) > 1 :
        repository = sys.argv[1]
    
    m_mu = 0.105 # GeV
    m_el = 0.000511 # GeV

    el_W = ( 'el_W_mA0p1_ebeam2p0.root' , 
        [('Scaled from 2.2 GeV', 'el_W_2p2to2.root'),
         ('Scaled from 2.5 GeV', 'el_W_2p5to2.root'),
         ('Scaled from 3.0 GeV', 'el_W_3to2.root'  ),
         ('Scaled from 4.0 GeV', 'el_W_4to2.root'  )
        ])

    el_Cu = ( 'el_Cu_mA0p1_ebeam2p0.root' , 
        [('Scaled from 2.2 GeV', 'el_Cu_2p2to2.root'),
         ('Scaled from 2.5 GeV', 'el_Cu_2p5to2.root'),
         ('Scaled from 3.0 GeV', 'el_Cu_3to2.root'  ),
         ('Scaled from 4.0 GeV', 'el_Cu_4to2.root'  )
        ])

    el_Si = ( 'el_Si_mA0p1_ebeam2p0.root' , 
        [('Scaled from 2.2 GeV', 'el_Si_2p2to2.root'),
         ('Scaled from 2.5 GeV', 'el_Si_2p5to2.root'),
         ('Scaled from 3.0 GeV', 'el_Si_3to2.root'  ),
         ('Scaled from 4.0 GeV', 'el_Si_4to2.root'  )
        ])

    mu_Pb = ( 'mu_Pb_mA1p0_ebeam100.root' ,
       [('Scaled from 111 GeV', 'mu_Pb_111to100.root'),
        ('Scaled from 125 GeV', 'mu_Pb_125to100.root'),
        ('Scaled from 150 GeV', 'mu_Pb_150to100.root'),
        ('Scaled from 200 GeV', 'mu_Pb_200to100.root')
       ])

    mu_Cu = ( 'mu_Cu_mA1p0_ebeam100.root' ,
       [('Scaled from 111 GeV', 'mu_Cu_111to100.root'),
        ('Scaled from 125 GeV', 'mu_Cu_125to100.root'),
        ('Scaled from 150 GeV', 'mu_Cu_150to100.root'),
        ('Scaled from 200 GeV', 'mu_Cu_200to100.root')
       ])

    
    scaling_validation(el_W, 'Electron' , 
                       'Electrons on Tungsten\n$m_{A\'} = 0.1$ GeV', 2., 
                       ke_ratio_ylim=(0.8,1.15), 
                       ang_ratio_ylim=(0.8,1.25), ang_raw_ylim=(4e-4,4.),
                       ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
                       file_prefix='el_W')
    
    scaling_validation(mu_Pb, 'Muon' , 
                       'Muons on Lead\n$m_{A\'} = 1.0$ GeV', 100.,
                       ke_ratio_ylim=(0.8,1.05),
                       ke_legend_kw=dict(bbox_to_anchor=(0.95,0.),loc='lower right'),
                       file_prefix='mu_Pb')
    
    efrac_pt_rates([
        ('Si Target, Z = 14', el_Si[0]),
        ('Cu Target, Z = 29', el_Cu[0]),
        ('W Target, Z = 74', el_W[0])],
        '2 GeV Electrons\n$m_{A\'} = 0.1$ GeV',
        2.0,
        file_prefix='el_2GeV')
    
if __name__ == '__main__' :
    main()
