"""Plot total xsec and their ratio to MadGraph"""

import pandas as pd
import matplotlib.pyplot as plt
import mplhep
plt.style.use(mplhep.style.ROOT)

import scipy

def xsec_plot(mg, others, title = None, ratio = True) :
    fig, (raw, ratio) = plt.subplots(ncols=1,nrows=2,sharex='col', 
                                     gridspec_kw = {'height_ratios' : [3,1]})
    plt.subplots_adjust(hspace=0.)

    raw.set_title(title)
    raw.plot(mg['Energy [GeV]'], mg['Xsec [pb]'], label='MG')
    ratio.plot(mg['Energy [GeV]'], [1. for i in range(len(mg['Energy [GeV]']))])
    for name, data in others :
        x = data['Energy [MeV]']/1000.
        y = data['Xsec [pb]']
        raw.plot(x, y, label=name)
        data_interp = scipy.interpolate.interp1d(x, y)
        data_at_mge = [data_interp(e) for e in mg['Energy [GeV]']]
        ratio.plot(mg['Energy [GeV]'], data_at_mge/mg['Xsec [pb]'])

    raw.set_ylabel('Total Cross Section [pb]')
    raw.legend()

    ratio.set_ylabel('Ratio to MG')
    ratio.set_xlabel('Incident Lepton Energy [GeV]')
    ratio.axhline((127/137)**3, color='gray')

def average_mg(all_samples) :
    """
    average the samples at each energy taking care to drop samples with less than 2/3 of the max

    this limit is just arbitrarily chosen to cut out outlier points that were seen 
    when plotting all samples as a scatter plot
    """

    return pd.read_csv(all_samples).groupby('Energy [GeV]').apply(lambda s : s[s > s.max()/1.5].mean())


def main() :
    xsec_plot(average_mg('data/mg/mu_xsec.csv'), [
            ('G4DB WW', pd.read_csv('data/dev/mu_xsec.csv'))
          ],
        title = 'Muons on Copper')
    plt.savefig('mu_xsec')

    xsec_plot(average_mg('data/mg/el_xsec.csv'), [
            ('G4DB WW', pd.read_csv('data/dev/el_xsec.csv')),
            ('G4DB IWW', pd.read_csv('data/dev/el_xsec_iww.csv'))
          ],
        title = 'Electrons on Tungsten')
    plt.savefig('el_xsec')

if __name__ == '__main__' :
    main()
