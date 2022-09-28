"""Analyze the dark brem events passed to us on the command line."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.ROOT)

import uproot
import numpy as np
import pandas as pd
import dark_brem_lhe
import scipy

mu_beam = 100105.658372 #MeV
el_beam = 4000.510999 #MeV

def read(beam_E, bias, fp) :
    """Read the passed file path into a dataframe one way or another
    
    Also calculate higher-level kinematic variables like recoil p_T,
    recoil lepton energy fraction, and recoil polar angle.
    
    If 'dblib' is in the file path, the `dark_brem_lhe` python module
    is used to parse the LHE file (if fp ends with '.lhe') or the directory
    of LHE files (otherwise) into the data frame.
    
    Otherwise, uproot is used to parse the ntuples at `dbint/dbint` into
    the data frame.
    """
    if 'dblib' in fp :
        if fp.endswith('.lhe') :
            df = dark_brem_lhe.DarkBremEventFile(fp).events
        else :
            df = dark_brem_lhe.DarkBremEventLibrary(fp).events()
        df['weight'] = np.ones(len(df['incident_mass']))/len(df['incident_mass'])
        df['incident_kinetic_energy_GeV'] = df['incident_energy']
        df['visible_energy'] = df['recoil_energy']
        df['beam_energy'] = beam_E
    else :
        with uproot.open(fp) as f :
            df = f['dbint/dbint'].arrays(library='pd')
        df['incident_kinetic_energy_GeV'] = (df['incident_energy']-df['incident_mass'])/1000.
    # some events in thin-target muon case have all kinematics in row set to DBL_MIN,
    #  this is due to when the framework "completes" an event without a successful simulation
    #  and therefore those rows can be dropped
    df.drop(df[df.weight < 1e-100].index, inplace=True)
    df['recoil_pt'] = (df['recoil_px']**2 + df['recoil_py']**2)**(1/2)
    df['energy_frac'] = (df['recoil_energy'] - df['incident_mass'])/(df['incident_energy']-df['incident_mass'])
    df['recoil_angle'] = np.arctan2(df['recoil_pt'],df['recoil_pz'])
    df['visible_energy_frac'] = df['visible_energy']/df['beam_energy']
    df['relative_weight'] = bias*df['weight']
    return df

def bundle(data_dir, mg_dir) :
    """Bundle the data from the input directory into stacked tuples+dicts to make plotting eaiser below
    
    This is where the names of the files are defined. If the output names are changed
    within sim/config.py, we would need to propagate that here. If the parameters that
    affect the output names are changed, we would need to propagate that here.
    """
    
    thin_tgt = (
        ('4 GeV Electrons\non 0.35 mm Tungsten',
         { # electrons
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_1.root'),
          'MG/ME' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
         }
        ),
        ('100 GeV Muons\non 100 mm Brass',
         { # muons
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_1.root'),
          'MG/ME' : read(mu_beam/1000.,5e5,f'{mg_dir}/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000_IncidentEnergy_100.0_unweighted_events.lhe')
         }
        )
    )
    
    thick_tgt = (
        ('4 GeV Electrons\non 18mm Tungsten',
         {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_1.root'),
          'Monoenergetic 4GeV MG/ME' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
         }
        ),
        ('100 GeV Muons\non 2m Brass',
         {
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_1.root'),
          'Monoenergetic 100GeV MG/ME' : read(mu_beam/1000.,5e5,f'{mg_dir}/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000_IncidentEnergy_100.0_unweighted_events.lhe')
         }
        )
    )

    na64 = ('100 GeV Electrons\non 1mm Lead',
        {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_lead_depthmm_1.0_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_lead_depthmm_1.0_mAMeV_100_events_50000_run_1.root')
        }
        )

    extra_thin = ('4 GeV Electrons\non 0.035mm Tungsten',
        {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_0.035_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_0.035_mAMeV_100_events_50000_run_1.root'),
          'MG' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
        }
        )
    
    return thin_tgt, thick_tgt, na64, extra_thin

def single(data_packet, kinematic_variable, xlabel, file_name,
           weight = True, ylabel = 'Weighted Event Fraction', yscale = 'log', 
           drop_mg = False, 
           hist_kwargs = {}, legend_kwargs = {}) :
    """Plot a single kinematic variable for the input data packet"""
    (title, data) = data_packet
    ((ax)) = plt.gcf().subplots()
    plt.gcf().set_size_inches(11,8)

    ax.set_xlabel(xlabel)
    ax.set_yscale(yscale)
    ax.set_ylabel(ylabel)
    for name, df in data.items() :
        if drop_mg and 'MG/ME' in name :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        ax.hist(df[kinematic_variable],
                weights = weights,
                label = name, linewidth = 2.,
                histtype = 'step', **hist_kwargs)
    l = ax.legend(title=title, **legend_kwargs)
    plt.setp(l.get_title(), multialignment='right')
    
    plt.savefig(file_name)
    plt.clf() 

def side_by_side(data_packet, kinematic_variable, xlabel, file_name,
                 weight = True, ylabel = 'Weighted Event Fraction', yscale = 'log', 
                 drop_mg = False,
                 el_kwargs = {}, mu_kwargs = {}, 
                 hist_kwargs = {}, legend_kwargs = {}) :
    """Plot same kinematic variable side-by-side with a shared y-axis"""
    
    ((el_title,el_data),(mu_title,mu_data)) = data_packet
    ((el_ax, mu_ax)) = plt.gcf().subplots(ncols = 2, nrows = 1, sharey = 'row')
    plt.gcf().set_size_inches(22,8)
    plt.subplots_adjust(hspace=0.3, wspace=0.)
    
    for ax in (el_ax,mu_ax) :
        ax.set_xlabel(xlabel)
        ax.set_yscale(yscale)
        
    el_ax.set_ylabel(ylabel)
    
    for name, df in el_data.items() :
        if drop_mg and 'MG/ME' in name :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        el_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **el_kwargs, **hist_kwargs)
    l = el_ax.legend(title=el_title+'\n$m_{A\'}=0.1$ GeV', **legend_kwargs)
    plt.setp(l.get_title(), multialignment='right')
    
    
    for name, df in mu_data.items() :
        if drop_mg and 'MG/ME' in name :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        mu_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **mu_kwargs, **hist_kwargs)
    l = mu_ax.legend(title = mu_title+'\n$m_{A\'}=1$ GeV', **legend_kwargs)
    plt.setp(l.get_title(), multialignment='right')
    
    plt.savefig(file_name)
    plt.clf()

def side_by_side_no_share(data_packet, kinematic_variable, xlabel, file_name,
                 weight = True, ylabel = 'Weighted Event Fraction', yscale = 'log', 
                 drop_mg = False, 
                 el_ylim = None, mu_ylim = None,
                 el_kwargs = {}, mu_kwargs = {}, 
                 hist_kwargs = {}, legend_kwargs = {}) :
    """Plot same kinematic variable side-by-side without a shared y-axis"""
    
    ((el_title,el_data),(mu_title,mu_data)) = data_packet
    ((el_ax, mu_ax)) = plt.gcf().subplots(ncols = 2, nrows = 1)
    plt.gcf().set_size_inches(22,8)
    
    for ax in (el_ax,mu_ax) :
        ax.set_xlabel(xlabel)
        ax.set_yscale(yscale)
        
    el_ax.set_ylabel(ylabel)
    
    for name, df in el_data.items() :
        if drop_mg and 'MG/ME' in name :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        el_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **el_kwargs, **hist_kwargs)
    if el_ylim is not None :
        el_ax.set_ylim(el_ylim)
    l = el_ax.legend(title=el_title+'\n$m_{A\'}=0.1$ GeV', **legend_kwargs)
    plt.setp(l.get_title(), multialignment='right')
    
    for name, df in mu_data.items() :
        if drop_mg and 'MG/ME' in name :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        mu_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **mu_kwargs, **hist_kwargs)
    if mu_ylim is not None :
        mu_ylim.set_ylim(mu_ylim)
    l = mu_ax.legend(title = mu_title+'\n$m_{A\'}=1$ GeV', **legend_kwargs)
    plt.setp(l.get_title(), multialignment='right')
    
    plt.savefig(file_name)
    plt.clf()

def xsec_plot(mg, others, file_name, title = None, ratio = True) :
    (raw, ratio) = plt.gcf().subplots(ncols = 1, nrows = 2, 
        sharex = 'col', gridspec_kw=dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace=0)

    mg_x = mg['Energy [GeV]']
    mg_y = mg['Xsec [pb]']*(127.9/137)**3

    raw.plot(mg_x, mg_y,marker='.', linewidth=0, label='MG')
    ratio.plot(mg_x, [1. for x in mg_x], marker='.', linewidth=0)
    for name, data in others :
        y = data['Xsec [pb]']
        if 'All MG' in name :
            x = data['Energy [GeV]']
            raw.plot(x,y*(127.9/137)**3,label=name,marker='.', linewidth=0)
        else :
            x = data['Energy [MeV]']/1000.
            raw.plot(x, y, label=name)
            data_interp = scipy.interpolate.interp1d(x, y)
            data_at_mge = [data_interp(e) for e in mg_x]
            ratio.plot(mg_x, data_at_mge/mg_y, marker='.', linewidth=0)

    raw.set_ylabel('Total Cross Section [pb]')
    l = raw.legend(title=title)
    plt.setp(l.get_title(), multialignment='right')

    ratio.set_ylabel('Ratio to MG')
    ratio.set_xlabel('Incident Lepton Energy [GeV]')
    plt.savefig(file_name)
    plt.clf()

def average_mg(all_samples) :
    """
    average the samples at each energy taking care to drop samples with less than 2/3 of the max

    this limit is just arbitrarily chosen to cut out outlier points that were seen 
    when plotting all samples as a scatter plot
    """

    return pd.read_csv(all_samples).groupby('Energy [GeV]').apply(lambda s : s[s > s.max()/1.5].mean())

        
def main() :
    import argparse
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir',help='Directory data is in')
    parser.add_argument('--out_dir',help='Directory to put plots (Default: data_dir)')
    parser.add_argument('--mg_dir',help='Directory MG libraries are in',default='dblib')
    parser.add_argument('--xsec-only', help='Only print xsec plots',action='store_true')
    
    arg = parser.parse_args()
    
    # make sure output directory exists
    if arg.out_dir is None :
        arg.out_dir = arg.data_dir
    os.makedirs(arg.out_dir, exist_ok=True)

    xsec_plot(average_mg('data/mg/mu_xsec.csv'), [
            ('DMG4 WW', pd.read_csv(f'{arg.data_dir}/dmg4_mu_xsec.csv')),
            ('G4DB WW', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec.csv')),
            ('All MG', pd.read_csv('data/mg/mu_xsec.csv'))
          ],
        f'{arg.out_dir}/mu_xsec.pdf',
        title = 'Muons on Copper')

    xsec_plot(average_mg('data/mg/el_xsec.csv'), [
            ('DMG4 IWW + Log Approx + K factors', 
              pd.read_csv(f'{arg.data_dir}/dmg4_el_xsec.csv').sort_values('Energy [MeV]')),
            #('G4DB WW', pd.read_csv('data/dev/el_xsec.csv')),
            ('G4DB IWW', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec.csv')),
            ('All MG', pd.read_csv('data/mg/el_xsec.csv'))
          ],
         f'{arg.out_dir}/el_xsec.pdf',
        title = 'Electrons on Tungsten')

    if arg.xsec_only :
        return

    # load data into memory bundles
    thin_tgt, thick_tgt, na64, extra_thin = bundle(arg.data_dir, arg.mg_dir)
    
    # get to plotting
    single(na64, 'recoil_angle', 'Lepton Recoil Angle [rad]',
           #el_ylim = (7e-4,2),
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = arg.out_dir+'/100GeV-electron-lead-recoil-angle.pdf')
    single(na64, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = arg.out_dir+'/100GeV-electron-lead-visible-energy.pdf')
    single(na64, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower center'},
           file_name = arg.out_dir+'/100GeV-electron-lead-visible-energy-cumulative.pdf')
    single(na64, 'incident_kinetic_energy_GeV', 
           'Lepton Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,100.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = arg.out_dir+'/100GeV-electron-lead-incident-energy.pdf')
    single(na64, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = arg.out_dir+'/100GeV-electron-lead-event-weight.pdf')

    single(extra_thin, 'recoil_angle', 'Lepton Recoil Angle [rad]',
           #el_ylim = (7e-4,2),
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = arg.out_dir+'/4GeV-electron-extra-thin-recoil-angle.pdf')
    single(extra_thin, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = arg.out_dir+'/4GeV-electron-extra-thin-visible-energy.pdf')
    single(extra_thin, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'framealpha': 0.8},
           file_name = arg.out_dir+'/4GeV-electron-extra-thin-visible-energy-cumulative.pdf')
    single(extra_thin, 'incident_kinetic_energy_GeV', 
           'Lepton Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,100.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = arg.out_dir+'/4GeV-electron-extra-thin-incident-energy.pdf')
    single(extra_thin, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = arg.out_dir+'/4GeV-electron-extra-thin-event-weight.pdf')

    side_by_side_no_share(thin_tgt, 'recoil_angle', 'Lepton Recoil Angle [rad]',
                 el_ylim = (7e-4,2),
                 hist_kwargs = {'range' : (0,2), 'bins' : 50},
                 file_name = arg.out_dir+'/thin-recoil-angle.pdf')
    side_by_side(thin_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
                 hist_kwargs = {'range' : (0,1), 'bins' : 50},
                 file_name = arg.out_dir+'/thin-visible-energy.pdf')
    side_by_side(thin_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
                 ylabel = 'Fraction Events Below Energy Cut',
                 yscale = 'linear',
                 hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
                 legend_kwargs = {'loc':'lower right', 'framealpha' : 0.8},
                 file_name = arg.out_dir+'/thin-visible-energy-cumulative.pdf')
    side_by_side(thin_tgt, 'incident_kinetic_energy_GeV', 
                 'Lepton Kinetic Energy Prior to DB [GeV]',
                 el_kwargs = {'range' : (0.,4.)},
                 mu_kwargs = {'range' : (0.,100.) },
                 legend_kwargs = {'loc' : 'upper left'},
                 hist_kwargs = { 'bins' : 50 },
                 drop_mg = True,
                 file_name = arg.out_dir+'/thin-incident-energy.pdf')
    side_by_side(thin_tgt, 'relative_weight', 'Event Weight',
                 weight = False, 
                 hist_kwargs = {'range':(1,1.2),'bins':50},
                 drop_mg = True,
                 file_name = arg.out_dir+'/thin-event-weight.pdf')
    
    side_by_side(thick_tgt, 'recoil_angle', 'Lepton Recoil Angle [rad]',
                 hist_kwargs = {'range' : (0,2), 'bins' : 50},
                 file_name = arg.out_dir+'/thick-recoil-angle.pdf')
    side_by_side(thick_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
                 hist_kwargs = {'range' : (0,1), 'bins' : 50},
                 file_name = arg.out_dir+'/thick-visible-energy.pdf')
    side_by_side(thick_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
                 ylabel = 'Fraction Events Below Energy Cut',
                 yscale = 'linear',
                 hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
                 legend_kwargs = {'loc':'lower center'},
                 file_name = arg.out_dir+'/thick-visible-energy-cumulative.pdf')
    side_by_side(thick_tgt, 'incident_kinetic_energy_GeV', 
                 'Lepton Kinetic Energy Prior to DB [GeV]',
                 el_kwargs = {'range' : (0.,4.)},
                 mu_kwargs = {'range' : (0.,100.) },
                 legend_kwargs = {'loc' : 'upper left'},
                 hist_kwargs = { 'bins' : 50 },
                 drop_mg = True,
                 file_name = arg.out_dir+'/thick-incident-energy.pdf')
    side_by_side(thick_tgt, 'relative_weight', 'Event Weight',
                 weight = False,
                 el_kwargs = {'range':(1,10), 'bins':100},
                 mu_kwargs = {'range':(1,1.2),'bins':50},
                 drop_mg = True,
                 file_name = arg.out_dir+'/thick-event-weight.pdf')
    
if __name__ == '__main__' :
    main()
