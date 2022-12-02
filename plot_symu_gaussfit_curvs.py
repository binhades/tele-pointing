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

def generate_simu_spec(fwhm=5,snr=3,nbin=20,seed=100):
    # snr -- signal to noise
    # Beam FWHM: 5 arcmin
    #fwhm = 5 # arcmin beam
    # scan length: 5 FWHM
    sleng = 5 * fwhm # 25 arcmin
    length = 5 * nbin

    x = (fwhm/nbin) * (np.arange(length) - length/2)
    np.random.seed(seed=seed)
    noi = np.random.normal(scale=1/snr,size=length)
    
    a = 1.0 # Jy
    b = 0 
    c = fwhm/(2.0 *(2.0*np.log(2))**0.5)
    g = Gaussian(x, a, b, c)
    
    power = g + noi
    
    return x, power, g, noi


def plot_gfit_curve(bins,snr,offset,offerr, fout=None,isshow=True,format='png',metadata=None):

    fig, axes = plt.subplots(1,2,sharex=True,figsize=(12,6),gridspec_kw={'height_ratios': [1]})

    axes[0].plot(snr, offset[0,:], color='dimgrey',  linestyle='--',     lw=2., label=f'N-Bins: {bins[0]:d}')
    axes[0].plot(snr, offset[1,:], color='k',        linestyle='dotted', lw=2., label=f'N-Bins: {bins[1]:d}')
    axes[0].plot(snr, offset[2,:], color='tab:red',  linestyle='-.',     lw=2., label=f'N-Bins: {bins[2]:d}')
    axes[0].plot(snr, offset[3,:], color='tab:blue', linestyle='-',      lw=2., label=f'N-Bins: {bins[3]:d}')

    axes[1].plot(snr, offerr[0,:], color='dimgrey',  linestyle='--',     lw=2., label=f'N-Bins: {bins[0]:d}')
    axes[1].plot(snr, offerr[1,:], color='k',        linestyle='dotted', lw=2., label=f'N-Bins: {bins[1]:d}')
    axes[1].plot(snr, offerr[2,:], color='tab:red',  linestyle='-.',     lw=2., label=f'N-Bins: {bins[2]:d}')
    axes[1].plot(snr, offerr[3,:], color='tab:blue', linestyle='-',      lw=2., label=f'N-Bins: {bins[3]:d}')
    #axes[0].plot(snr, y_gauss, color='tab:blue',  linestyle='-',     lw=2.0, label='Gaussian Fit')
    #axes[1].plot(snr, y_nois,  color='grey',     linestyle='--',      lw=1.0, label='Residual') 

    h0,l0 = axes[0].get_legend_handles_labels()
    h1,l1 = axes[1].get_legend_handles_labels()
    axes[0].legend(h0,l0, loc='upper right',fontsize=14) 
    axes[1].legend(h1,l1, loc='upper right',fontsize=14)

    axes[0].tick_params(labelsize=18)
    axes[1].tick_params(labelsize=18)

    axes[0].set_xlabel('Signal-to-Noise Ratio',fontsize=18)
    axes[1].set_xlabel('Signal-to-Noise Ratio',fontsize=18)
    axes[0].set_ylabel('Fitting Offsets (arcmin)',fontsize=18)
    axes[1].set_ylabel('Fitting Errors (arcmin)',fontsize=18)

    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    if fout is not None:
        plt.savefig(fout,format=format,dpi=300,bbox_inches='tight',metadata=metadata)
    if isshow:
        plt.show()

    plt.close(fig)

    return 0

def main(args):
    offset = []
    offerr = []
    bin_list = [10,20,40,80]
    snr_list = np.arange(3,50)
    for nbin in bin_list:
        offset_l = []
        offerr_l = []
        for snr in snr_list:
            #--------------------------------------------------------
            # Simulat signal
            x, y, y_gauss, y_nois = generate_simu_spec(fwhm=args.fwhm, snr=snr,nbin=nbin)
            #--------------------------------------------------------
            # Gaussian Fit
            fitfunc,pfit,efit = fit_gauss(x,y)
            #--------------------------------------------------------
            offset_l.append(abs(pfit[1]))
            offerr_l.append(efit[1,1])
        offset.append(offset_l)
        offerr.append(offerr_l)

    offset = np.array(offset)
    offerr = np.array(offerr)
    print(offset.shape)
    print(snr_list.shape)

    #--------------------------------------------------------
    # Plots
    plot_gfit_curve(bin_list, snr_list, offset, offerr, isshow=args.isshow,fout=args.fileout)
    #--------------------------------------------------------

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fileout', type=str, help='the output file name')
    parser.add_argument('--isshow',  action='store_true', help='the output file name')
    parser.add_argument('--fwhm',    type=float, default=5, help='beam size in arcmin, 5')
    args = parser.parse_args()
    main(args)
