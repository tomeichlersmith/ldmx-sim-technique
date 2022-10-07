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
    
    thin_el = ('thin-electron',
          '$m_{A\'} = 0.1$ GeV\n4 GeV Electrons\non 0.35 mm Tungsten',
         { # electrons
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_0.35_mAMeV_100_events_50000_run_1.root'),
          'MG/ME' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
         }
        )

    thin_mu = ('thin-muon',
          '$m_{A\'} = 1$ GeV\n100 GeV Muons\non 100 mm Brass',
         { # muons
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_100.0_mAMeV_1000_events_50000_run_1.root'),
          'MG/ME' : read(mu_beam/1000.,5e5,f'{mg_dir}/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000_IncidentEnergy_100.0_unweighted_events.lhe')
         }
        )
    
    thick_el = ('thick-electron',
          '$m_{A\'} = 0.1$ GeV\n4 GeV Electrons\non 18mm Tungsten',
         {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_18.0_mAMeV_100_events_50000_run_1.root'),
          'Monoenergetic 4GeV MG/ME' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
         }
        )

    thick_mu = ('thick-muon',
          '$m_{A\'} = 1$ GeV\n100 GeV Muons\non 2m Brass',
         {
          'G4DarkBreM' : read(mu_beam,1e7,f'{data_dir}/ntuple_g4db_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_3000.root'),
          'DMG4' : read(mu_beam,1e11,f'{data_dir}/ntuple_dmg4_muon_brass_depthmm_2000.0_mAMeV_1000_events_50000_run_1.root'),
          'Monoenergetic 100GeV MG/ME' : read(mu_beam/1000.,5e5,f'{mg_dir}/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000/muon_copper_MaxE_100.0_MinE_2.0_RelEStep_0.1_UndecayedAP_mA_1.0_run_3000_IncidentEnergy_100.0_unweighted_events.lhe')
         }
        )

    na64 = ('100GeV-electron-lead',
        '$m_{A\'} = 0.1$ GeV\n100 GeV Electrons\non 1mm Lead',
        {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_lead_depthmm_1.0_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_lead_depthmm_1.0_mAMeV_100_events_50000_run_1.root')
        }
        )

    extra_thin = ('4GeV-electron-extra-thin',
        '$m_{A\'} = 0.1$ GeV\n4 GeV Electrons\non 0.035mm Tungsten',
        {
          'G4DarkBreM' : read(el_beam,1e8,f'{data_dir}/ntuple_g4db_electron_tungsten_depthmm_0.035_mAMeV_100_events_50000_run_3000.root'),
          'DMG4' : read(el_beam,1e12,f'{data_dir}/ntuple_dmg4_electron_tungsten_depthmm_0.035_mAMeV_100_events_50000_run_1.root'),
          'MG' : read(el_beam/1000.,5e5,f'{mg_dir}/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000/electron_tungsten_MaxE_4.0_MinE_0.2_RelEStep_0.1_UndecayedAP_mA_0.1_run_3000_IncidentEnergy_4.0_unweighted_events.lhe')
        }
        )
    
    return thin_el, thin_mu, thick_el, thick_mu, na64, extra_thin

def single(data_packet, kinematic_variable, xlabel, file_name,
           weight = True, ylabel = 'Weighted Event Fraction', yscale = 'log', 
           drop_mg = False, ylim = None,
           hist_kwargs = {}, legend_kwargs = {}) :
    """Plot a single kinematic variable for the input data packet"""
    (_, title, data) = data_packet
    ((ax)) = plt.gcf().subplots()
    plt.gcf().set_size_inches(11,8)

    ax.set_xlabel(xlabel)
    ax.set_yscale(yscale)
    ax.set_ylabel(ylabel)
    if ylim is not None :
        ax.set_ylim(ylim)
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

def xsec_plot(mg_csv, others, file_name, 
    xlabel = 'Incident Lepton Energy [GeV]', title = None) :
    (raw, ratio) = plt.gcf().subplots(ncols = 1, nrows = 2, 
        sharex = 'col', gridspec_kw=dict(height_ratios = [3,1]))
    plt.subplots_adjust(hspace=0)

    mg = pd.read_csv(mg_csv).groupby('Energy [GeV]').apply(lambda samples : samples[samples > samples.median() - 2*samples.std()].mean()).drop(columns='Energy [GeV]').reset_index()

    mg_x = mg['Energy [GeV]']
    mg_y = mg['Xsec [pb]']*(127.9/137)**3

    for name, data in others :
        y = data['Xsec [pb]']
        x = data['Energy [MeV]']/1000.
        raw.plot(x, y, label=name)
        data_interp = scipy.interpolate.interp1d(x, y)
        data_at_mge = [data_interp(e) for e in mg_x]
        ratio.plot(mg_x, data_at_mge/mg_y, marker='.', linewidth=0)

    raw.plot(mg_x, mg_y,marker='.', linewidth=0, label='MG/ME')
    ratio.plot(mg_x, [1. for x in mg_x], marker='.', markersize=0, linewidth=0)

    raw.set_ylabel('Total Cross Section / $\epsilon^2$ [pb]')
    l = raw.legend(title=title)
    plt.setp(l.get_title(), multialignment='right')

    ratio.set_ylabel('Ratio to MG/ME')
    ratio.set_xlabel(xlabel)
    plt.savefig(file_name)
    plt.clf()
        
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

    xsec_plot('data/mg/mu_xsec.csv', [
            ('G4DarkBreM', pd.read_csv(f'{arg.data_dir}/g4db_mu_xsec.csv')),
          ],
        f'{arg.out_dir}/mu_xsec.pdf',
        xlabel = 'Incident Muon Energy [GeV]',
        title = '$m_{A\'} = 0.2$ GeV\nMuons on Copper')

    xsec_plot('data/mg/el_xsec.csv', [
            ('G4DarkBreM', pd.read_csv(f'{arg.data_dir}/g4db_el_xsec.csv')),
          ],
         f'{arg.out_dir}/el_xsec.pdf',
        xlabel = 'Incident Electron Energy [GeV]',
        title = '$m_{A\'} = 0.1$ GeV\nElectrons on Tungsten')

    if arg.xsec_only :
        return

    # load data into memory bundles
    thin_el, thin_mu, thick_el, thick_mu, na64, extra_thin = bundle(arg.data_dir, arg.mg_dir)

    def filename(prefix, tail) :
        return f'{arg.out_dir}/{prefix}-{tail}.pdf'
    
    # get to plotting
    single(na64, 'recoil_angle', 'Electron Recoil Angle [rad]',
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(na64[0],'recoil-angle'))
    single(na64, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(na64[0],'visible-energy'))
    single(na64, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower center'},
           file_name = filename(na64[0],'visible-energy-cumulative'))
    single(na64, 'incident_kinetic_energy_GeV', 
           'Electron Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,100.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(na64[0],'incident-energy'))
    single(na64, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = filename(na64[0],'event-weight'))

    single(extra_thin, 'recoil_angle', 'Lepton Recoil Angle [rad]',
           #el_ylim = (7e-4,2),
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(extra_thin[0],'recoil-angle'))
    single(extra_thin, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(extra_thin[0], 'visible-energy'))
    single(extra_thin, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'framealpha': 0.8},
           file_name = filename(extra_thin[0], 'visible-energy-cumulative'))
    single(extra_thin, 'incident_kinetic_energy_GeV', 
           'Lepton Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,4.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(extra_thin[0], 'incident-energy'))
    single(extra_thin, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = filename(extra_thin[0], 'event-weight'))

    single(thin_el, 'recoil_angle', 'Electron Recoil Angle [rad]',
           ylim = (7e-4,2),
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(thin_el[0],'recoil-angle'))
    single(thin_el, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(thin_el[0], 'visible-energy'))
    single(thin_el, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'bbox_to_anchor':(0.95,0.)},
           file_name = filename(thin_el[0], 'visible-energy-cumulative'))
    single(thin_el, 'incident_kinetic_energy_GeV', 
           'Electron Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,4.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(thin_el[0], 'incident-energy'))
    single(thin_el, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = filename(thin_el[0], 'event-weight'))

    single(thin_mu, 'recoil_angle', 'Muon Recoil Angle [rad]',
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(thin_mu[0],'recoil-angle'))
    single(thin_mu, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(thin_mu[0], 'visible-energy'))
    single(thin_mu, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'bbox_to_anchor':(0.95,0)},
           file_name = filename(thin_mu[0], 'visible-energy-cumulative'))
    single(thin_mu, 'incident_kinetic_energy_GeV', 
           'Muon Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,100.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(thin_mu[0], 'incident-energy'))
    single(thin_mu, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = filename(thin_mu[0], 'event-weight'))

    single(thick_el, 'recoil_angle', 'Electron Recoil Angle [rad]',
           ylim = (7e-4,2),
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(thick_el[0],'recoil-angle'))
    single(thick_el, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(thick_el[0], 'visible-energy'))
    single(thick_el, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'bbox_to_anchor':(0.95,0)},
           file_name = filename(thick_el[0], 'visible-energy-cumulative'))
    single(thick_el, 'incident_kinetic_energy_GeV', 
           'Electron Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,4.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(thick_el[0], 'incident-energy'))
    single(thick_el, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,10.),'bins':50},
           drop_mg = True,
           file_name = filename(thick_el[0], 'event-weight'))

    single(thick_mu, 'recoil_angle', 'Muon Recoil Angle [rad]',
           hist_kwargs = {'range' : (0,2), 'bins' : 50},
           file_name = filename(thick_mu[0],'recoil-angle'))
    single(thick_mu, 'visible_energy_frac', 'Visible Energy Fraction of Beam',
           hist_kwargs = {'range' : (0,1), 'bins' : 50},
           file_name = filename(thick_mu[0], 'visible-energy'))
    single(thick_mu, 'visible_energy_frac', 'Visible Energy Fraction of Beam', 
           ylabel = 'Fraction Events Below Energy Cut',
           yscale = 'linear',
           hist_kwargs = {'range' : (0,1), 'bins': 50, 'cumulative' : True},
           legend_kwargs = {'loc':'lower right', 'bbox_to_anchor':(0.95,0)},
           file_name = filename(thick_mu[0], 'visible-energy-cumulative'))
    single(thick_mu, 'incident_kinetic_energy_GeV', 
           'Muon Kinetic Energy Prior to DB [GeV]',
           hist_kwargs = {'range' : (0.,100.), 'bins' : 50 },
           legend_kwargs = {'loc' : 'upper left'},
           drop_mg = True,
           file_name = filename(thick_mu[0], 'incident-energy'))
    single(thick_mu, 'relative_weight', 'Event Weight',
           weight = False, 
           hist_kwargs = {'range':(1,1.2),'bins':50},
           drop_mg = True,
           file_name = filename(thick_mu[0], 'event-weight'))

if __name__ == '__main__' :
    main()
