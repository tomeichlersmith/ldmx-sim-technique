"""Analyze the dark brem events passed to us on the command line."""

import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot
import mplhep
import numpy as np
import dark_brem_lhe

plt.style.use(mplhep.style.ROOT)
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
        ('4 GeV Electrons on 0.35 mm Tungsten',
         { # electrons
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_1.root'),
          'MG' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
         }
        ),
        ('100 GeV Muons on 100 mm Brass',
         { # muons
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_1.root'),
          'MG' : read(mu_beam/1000.,5e5,f'{mg_dir}/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000_IncidentEnergy_100.0_unweighted_events.lhe')
         }
        )
    )
    
    thick_tgt = (
        ('4 GeV Electrons on 18mm Tungsten',
         {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_1.root')
         }
        ),
        ('100 GeV Muons on 2m Brass',
         {
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_1.root')
         }
        )
    )
    
    return thin_tgt, thick_tgt

def side_by_side(data_packet, kinematic_variable, xlabel, file_name,
                 weight = True, ylabel = 'Weighted Event Fraction', yscale = 'log', 
                 drop_mg = False, fig_ext = ['pdf'],
                 el_kwargs = {}, mu_kwargs = {}, 
                 hist_kwargs = {}, legend_kwargs = {}) :
    """Plot same kinematic variable side-by-side with a shared y-axis"""
    
    ((el_title,el_data),(mu_title,mu_data)) = data_packet
    fig, ((el_ax, mu_ax)) = plt.subplots(ncols = 2, nrows = 1, sharey = 'row')
    fig.set_size_inches(22,8)
    plt.subplots_adjust(hspace=0.3, wspace=0.)
    
    for ax in (el_ax,mu_ax) :
        ax.set_xlabel(xlabel)
        ax.set_yscale(yscale)
        
    el_ax.set_ylabel(ylabel)
    
    el_ax.set_title(el_title)
    for name, df in el_data.items() :
        if drop_mg and name == 'MG' :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        el_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **el_kwargs, **hist_kwargs)
    el_ax.legend(title='$m_{A\'}=0.1$ GeV', **legend_kwargs)
    
    mu_ax.set_title(mu_title)
    for name, df in mu_data.items() :
        if drop_mg and name == 'MG' :
            continue
        weights = None
        if weight :
            weights = df['weight']/df['weight'].sum()
        mu_ax.hist(df[kinematic_variable],
                   weights = weights,
                   label = name, linewidth = 2.,
                   histtype = 'step', **mu_kwargs, **hist_kwargs)
    mu_ax.legend(title = '$m_{A\'}=1$ GeV', **legend_kwargs)
    
    plt.savefig(file_name)
        
def main() :
    import argparse
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument('data_dir',help='Directory data is in')
    parser.add_argument('--out_dir',help='Directory to put plots (Default: data_dir)')
    parser.add_argument('--mg_dir',help='Directory MG libraries are in',default='dblib')
    
    arg = parser.parse_args()
    
    # load data into memory bundles
    thin_tgt, thick_tgt = bundle(arg.data_dir, arg.mg_dir)
    
    # make sure output directory exists
    if arg.out_dir is None :
        arg.out_dir = arg.data_dir
    os.makedirs(arg.out_dir, exist_ok=True)
    
    # get to plotting
    side_by_side(thin_tgt, 'recoil_angle', 'Lepton Recoil Angle [rad]',
                 hist_kwargs = {'range' : (0,2), 'bins' : 50},
                 file_name = arg.out_dir+'/thin-recoil-angle')
    side_by_side(thin_tgt, 'visible_energy_frac', 'Visible Energy Fraciton of Beam',
                 hist_kwargs = {'range' : (0,1), 'bins' : 50},
                 file_name = arg.out_dir+'/thin-visible-energy')
    side_by_side(thin_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
                 ylabel = 'Fraction Events Below Energy Cut',
                 hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
                 legend_kwargs = {'loc':'lower center'},
                 file_name = arg.out_dir+'/thin-visible-energy-cumulative')
    side_by_side(thin_tgt, 'incident_kinetic_energy_GeV', 
                 'Lepton Kinetic Energy Prior to DB [GeV]',
                 el_kwargs = {'range' : (0.,4.)},
                 mu_kwargs = {'range' : (0.,100.) },
                 legend_kwargs = {'loc' : 'upper left'},
                 hist_kwargs = { 'bins' : 50 },
                 drop_mg = True,
                 file_name = arg.out_dir+'/thin-incident-energy')
    side_by_side(thin_tgt, 'relative_weight', 'Event Weight',
                 weight = False, 
                 hist_kwargs = {'range':(1,1.2),'bins':50},
                 drop_mg = True,
                 file_name = arg.out_dir+'/thin-event-weight')
    
    side_by_side(thick_tgt, 'recoil_angle', 'Lepton Recoil Angle [rad]',
                 hist_kwargs = {'range' : (0,2), 'bins' : 50},
                 file_name = arg.out_dir+'/thick-recoil-angle')
    side_by_side(thick_tgt, 'visible_energy_frac', 'Visible Energy Fraciton of Beam',
                 hist_kwargs = {'range' : (0,1), 'bins' : 50},
                 file_name = arg.out_dir+'/thick-visible-energy')
    side_by_side(thick_tgt, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
                 ylabel = 'Fraction Events Below Energy Cut',
                 hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
                 legend_kwargs = {'loc':'lower center'},
                 file_name = arg.out_dir+'/thick-visible-energy-cumulative')
    side_by_side(thick_tgt, 'incident_kinetic_energy_GeV', 
                 'Lepton Kinetic Energy Prior to DB [GeV]',
                 el_kwargs = {'range' : (0.,4.)},
                 mu_kwargs = {'range' : (0.,100.) },
                 legend_kwargs = {'loc' : 'upper left'},
                 hist_kwargs = { 'bins' : 50 },
                 drop_mg = True,
                 file_name = arg.out_dir+'/thick-incident-energy')
    side_by_side(thick_tgt, 'relative_weight', 'Event Weight',
                 weight = False,
                 el_kwargs = {'range':(1,10), 'bins':100},
                 mu_kwargs = {'range':(1,1.2),'bins':50},
                 drop_mg = True,
                 file_name = arg.out_dir+'/thick-event-weight')
    
if __name__ == '__main__' :
    main()
