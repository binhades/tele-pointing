#!//usr/bin/python3

import numpy as np
import argparse
from scipy import optimize
import matplotlib.pyplot as plt

def Gaussian(x, a, b, c):
    return a * np.exp(-0.5*(x-b)**2/c**2)

def fit_gauss(xdata,ydata):
# FWHM = 2*(2ln2)^0.5*sigma=2.35482*sigma

    nbins = len(xdata)
    imax = np.argmax(ydata)
    pinit = np.array([ydata[imax],0,5])

    def fitfunc(p,x):
        retval = Gaussian(x,p[0],p[1],p[2])
        return retval
    def errfunc(p,x,y):
        return y - fitfunc(p,x)

    out = optimize.leastsq(errfunc,pinit,args=(xdata,ydata),full_output=True)
    pfit = out[0]
    efit = out[1]
    return fitfunc,pfit,efit

def plot(x, y, y_gauss, y_nois, fout=None,isshow=True,format='png', metadata=None):

    fig, axes = plt.subplots(2,1,sharex=True,figsize=(12,8),gridspec_kw={'height_ratios': [2, 1]})

    axes[0].step(x, y,       color='dimgrey',  linestyle='dotted', lw=2., label='Observation')
    axes[0].plot(x, y_gauss, color='tab:blue',  linestyle='-',     lw=2.0, label='Gaussian Fit')
    axes[1].plot(x, y_nois,  color='grey',     linestyle='--',      lw=1.0, label='Residual') 

    h0,l0 = axes[0].get_legend_handles_labels()
    h1,l1 = axes[1].get_legend_handles_labels()
    axes[0].legend(h0,l0, loc='upper right',fontsize=14) 
    #axes[1].legend(h1,l1, loc='upper right',fontsize=14)

    axes[0].tick_params(labelsize=18)
    axes[1].tick_params(labelsize=18)

    axes[0].set_ylabel('Relative intensity',fontsize=18)
    axes[1].set_ylabel('Residual',fontsize=18)
    axes[1].set_xlabel('Scan Offset (arcmin)',fontsize=18)

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    if fout is not None:
        plt.savefig(fout,format=format,dpi=300,bbox_inches='tight',metadata=metadata)
    if isshow:
        plt.show()

    plt.close(fig)

    return 0

def generate_simu_spec(fwhm=5,snr=3,nbin=20):
    # snr -- signal to noise
    # Beam FWHM: 5 arcmin
    #fwhm = 5 # arcmin beam
    # scan length: 5 FWHM
    sleng = 5 * fwhm # 25 arcmin
    length = 5 * nbin

    x = (fwhm/nbin) * (np.arange(length) - length/2)
    np.random.seed(seed=100)
    noi = np.random.normal(scale=1/snr,size=length)
    
    a = 1.0 # Jy
    b = 0 
    c = fwhm/(2.0 *(2.0*np.log(2))**0.5)
    g = Gaussian(x, a, b, c)
    
    power = g + noi
    
    return x, power, g, noi

def main(args):
    #--------------------------------------------------------
    # Simulat signal
    # scan speed
    spd = args.spd/60.0 # arcsec to arcmin [per sec]
    if args.nbin is None:
        nbin = int(args.fwhm/(spd*args.tau))
    else:
        nbin = args.nbin

    x, y, y_gauss, y_nois = generate_simu_spec(fwhm=args.fwhm, snr=args.snr,nbin=nbin)
    #--------------------------------------------------------
    mStr = 'Signal: fwhm={fwhm:.1f}, snr={snr:d}, nbin={nbin:d}'.format(
            fwhm=args.fwhm, snr=args.snr,nbin=nbin)
    print(mStr)
    #--------------------------------------------------------
    fitfunc,pfit,efit = fit_gauss(x,y)
    y_gfit = fitfunc(pfit,x)
    y_resi = y - y_gfit
    rms = np.std(y_resi)
    print(pfit)
    
    pfwhm = pfit[2] * (2.0 *(2.0*np.log(2))**0.5)
    efwhm = efit[2,2] * (2.0 *(2.0*np.log(2))**0.5)
    print(f'Gaussi: Peak; Offset; FWHM;')
    print(f'FitPar: {pfit[0]:.2f}; {pfit[1]:.2f};   {pfwhm:.2f};')
    print(f'Errors: {efit[0,0]:.2f}; {efit[1,1]:.2f};   {efwhm:.2f};')
    print(f'RMS: {rms:.2f}')

    #--------------------------------------------------------
    # Plots
    plot(x, y, y_gfit, y_resi, isshow=args.isshow,fout=args.fileout,metadata={'Title':mStr})
    #--------------------------------------------------------

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fileout', type=str, help='the output file name')
    parser.add_argument('--isshow',  action='store_true', help='the output file name')
    parser.add_argument('--fwhm',    type=float, default=5, help='beam size in arcmin, 5')
    parser.add_argument('--spd',     type=float, default=15, help='the scan speed in arcsec, 15')
    parser.add_argument('--tau',     type=float, default=1, help='the sample integral time, 1 s')
    parser.add_argument('--snr',     type=int,   default=5, help='the snr of simulated signal')
    parser.add_argument('--nbin',     type=int,  help='the snr of simulated signal')
    args = parser.parse_args()
    main(args)
