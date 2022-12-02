#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
===================================================================
To extract sources and convert to csv from VLA calibrator list
===================================================================
bla bla ~~~

*By: Bin Liu*

*License: BSD*

"""
import argparse, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import astropy.units as u

def load_catalog(file_cata):
    with open(file_cata,'r') as fin:
        hdr = fin.readline().rstrip()
    hdr = hdr.split(',')
    sour_list = np.genfromtxt(file_cata,delimiter=',',skip_header=1,dtype=None,names=hdr,encoding=None)
    return sour_list

def plot_hist(flux,n_bins=200, file_out=None,fext='png',isshow=True,plottype='line',xmin=0,xmax=5,fontsize='xx-large',xlabel='Flux [Jy]',ylabel='Counts'):

    fig, axs = plt.subplots(1, 2, tight_layout=True, figsize=(12,5))
    
    # N is the count in each bin, bins is the lower-limit of the bin
    N, bins, patches = axs[0].hist(flux, bins=n_bins)
    
    # We'll color code by height, but you could use any scalar
    fracs = N / N.max()
    
    # we need to normalize the data to 0..1 for the full range of the colormap
    norm = colors.Normalize(fracs.min(), fracs.max())
    
    # Now, we'll loop through our objects and set the color of each accordingly
    for thisfrac, thispatch in zip(fracs, patches):
        color = plt.cm.viridis(norm(thisfrac))
        thispatch.set_facecolor(color)
    
    # We can also normalize our inputs by the total number of counts
    axs[1].hist(flux, bins=n_bins, density=True)
    
    # Now we format the y-axis to display percentage
    axs[1].yaxis.set_major_formatter(PercentFormatter(xmax=5))
    axs[0].set_xlim([xmin,xmax])
    axs[1].set_xlim([xmin,xmax])
    axs[0].set_xlabel(xlabel,fontsize=fontsize)
    axs[1].set_xlabel(xlabel,fontsize=fontsize)
    axs[0].set_ylabel(ylabel,fontsize=fontsize)
    axs[0].tick_params(axis='both', labelsize=fontsize) #which='minor'
    axs[1].tick_params(axis='both', labelsize=fontsize) #which='minor'
    fig.suptitle('VLA Ku-band (2cm) Calibrator Flux Distribution',fontsize=fontsize)
    if file_out is not None:
        plt.savefig(file_out,format=fext,dpi=300)
    if isshow:
        plt.show()
    plt.close()
    return 0
#===================
def plot_equ_sky_scatter(sources,size=None,colorInd=None,file_out=None,isshow=True,fontsize='x-large',plotGP=True,ax_factor=1):

    fig = plt.figure(figsize=(12,7))
    ax = fig.add_subplot(111, projection="aitoff")

    if plotGP:
        l = np.arange(3600)/10.0
        b0 = np.zeros(3600)
        bp = +10 * np.ones(3600)
        bm = -10 * np.ones(3600)
        gc = SkyCoord(0,0, frame='galactic', unit=u.deg)
        gc_equ = gc.transform_to('icrs')
        
        gal0 = SkyCoord(l[:], b0[:], frame='galactic', unit=u.deg)
        galp = SkyCoord(l[:], bp[:], frame='galactic', unit=u.deg)
        galm = SkyCoord(l[:], bm[:], frame='galactic', unit=u.deg)
        equ0 = gal0.icrs
        equp = galp.transform_to('icrs')
        equm = galm.transform_to('icrs')
        ax.scatter(ax_factor*gc_equ.ra.wrap_at(180*u.deg).radian, gc_equ.dec.wrap_at(180*u.deg).radian,marker='*',s=100,color='black')
        ax.plot(ax_factor*equ0.ra.wrap_at(180*u.deg).radian, equ0.dec.wrap_at(180*u.deg).radian, alpha=0.5,color='green',linestyle='-')
        ax.plot(ax_factor*equp.ra.wrap_at(180*u.deg).radian, equp.dec.wrap_at(180*u.deg).radian, alpha=0.5,color='tab:green',linestyle='--')
        ax.plot(ax_factor*equm.ra.wrap_at(180*u.deg).radian, equm.dec.wrap_at(180*u.deg).radian, alpha=0.5,color='tab:green',linestyle='--')
       

#    ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian, sources.dec.wrap_at(180*u.deg).radian,s=size*2)
    ind0 = np.where(size < 0.5)
    ind1 = np.where((size >= 0.5) & (size < 1.0))
    ind2 = np.where((size >= 1.0) & (size < 3.0))
    ind3 = np.where(size >= 3.0)
    print(f'Number of Weak: {ind0[0].shape[0]}')
    print(f'Number of Middle: {ind1[0].shape[0]}')
    print(f'Number of Strong: {ind2[0].shape[0]}')
    print(f'Number of Insane: {ind3[0].shape[0]}')
    ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian[ind0], sources.dec.wrap_at(180*u.deg).radian[ind0],s=size[ind0]*20,marker='o', color='blue')
    ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian[ind1], sources.dec.wrap_at(180*u.deg).radian[ind1],s=size[ind1]*20,marker='<', color='red')
    ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian[ind2], sources.dec.wrap_at(180*u.deg).radian[ind2],s=size[ind2]*20,marker='s', color='tab:brown')
    ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian[ind3], sources.dec.wrap_at(180*u.deg).radian[ind3],s=size[ind3]*20,marker='x', color='purple')
    # -------------------------------------------------------------------------
    xt0 = ax_factor*np.array([-135,-90,-45,0,45,90,135])*np.pi/180
    xtl0 =         ['15$^h$','18$^h$','21$^h$','0$^h$','3$^h$','6$^h$','9$^h$']
    xt1 = ax_factor*np.array([-150,-120,-90,-60,-30,0,30,60,90,120,150])*np.pi/180
    xtl1 =         ['14$^h$','16$^h$','18$^h$','20$^h$','22$^h$','0$^h$','2$^h$','4$^h$','6$^h$','8$^h$','10$^h$']
    yt = np.array([-75,-60,-45,-30,-15,0,15,30,45,60,75])*np.pi/180
    ytl = ['-75$^{\circ}$','-60$^{\circ}$','-45$^{\circ}$','-30$^{\circ}$','-15$^{\circ}$','0$^{\circ}$','15$^{\circ}$','30$^{\circ}$','45$^{\circ}$','60$^{\circ}$','']
    ax.set_xticks(xt1)
    ax.set_xticklabels(xtl1,fontsize=fontsize)
    ax.set_yticks(yt)
    ax.set_yticklabels(ytl,fontsize=fontsize)
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.grid(color='gray', linestyle='dotted')
    ax.set_title('VLA Ku-band Calibrators',fontsize=fontsize)
    ax.set_xlabel('RA (J2000)',fontsize=fontsize)
    ax.set_ylabel('Dec (J2000)',fontsize=fontsize)
#    ax.set_facecolor('xkcd:battleship grey')
    # -------------------------------------------------------------------------
    plt.tight_layout()
    fig.patch.set_facecolor('white')
    if file_out is not None:
        fig.savefig(file_out, dpi=300,format='png',bbox_inches='tight')
    if isshow:
        plt.show()
    return 0
#===================
def catalog_to_Skycoord(sour_list):
    sour_coor = SkyCoord(ra=sour_list['RA'],dec=sour_list['DEC'],frame='icrs')
    return sour_coor
#===================
def flux_to_steps(flux,steps=[]):
    indx = np.empty(flux.shape[0])
    for i, s in enumerate(steps): # small -> large
        if i == 0:
            indx[np.where(flux<s)] = i
        indx[np.where(flux>=s)] = i + 1
    return indx.astype(int)
#===================
def main(args):
    # check if file exist
    if not os.path.isfile(args.file_csv):
        print(f"File {args.file_csv} does not exist.")
        return 0

    sour_list = load_catalog(args.file_csv)
    sour_coor = catalog_to_Skycoord(sour_list)
    indx = flux_to_steps(sour_list['Flux'],steps=[0.5,1])
    plot_equ_sky_scatter(sour_coor,size=sour_list['Flux'],colorInd=indx,\
            file_out=args.file_out,fontsize='xx-large', \
            ax_factor=args.ax_factor, isshow=args.isshow, plotGP=args.plotGP)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
            description='''To plot flux histgram''',
            epilog='version 0.1.')

    parser.add_argument('file_csv', type=str, help='the input catalog csv file') 
    parser.add_argument('--file_out', type=str, help='the output figure name')
    parser.add_argument('--ax_factor', type=int, default=-1,help='the output figure name')
    parser.add_argument('--plotGP', action='store_true', help='set to plot Galactic plane')
    parser.add_argument('--isshow', action='store_true', help='set to plot Galactic plane')

    args = parser.parse_args()
    main(args)

