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
from astropy.coordinates import SkyCoord
import astropy.units as u

#===================
def load_catalog_np(file_cata):
    sour_list = np.genfromtxt(file_cata,delimiter=',',dtype=None,names=True,encoding=None)
    return sour_list

#===================
def plot_equ_sky_scatter(sources,size=None,colorInd=None,file_out=None,isshow=True,fontsize='x-large',plotGP=True,ax_factor=-1,groups=None):

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
        ax.plot(ax_factor*equ0.ra.wrap_at(180*u.deg).radian, equ0.dec.wrap_at(180*u.deg).radian, alpha=0.5, zorder=1, color='green',     linestyle='-' )
        ax.plot(ax_factor*equp.ra.wrap_at(180*u.deg).radian, equp.dec.wrap_at(180*u.deg).radian, alpha=0.5, zorder=1, color='tab:green', linestyle='--')
        ax.plot(ax_factor*equm.ra.wrap_at(180*u.deg).radian, equm.dec.wrap_at(180*u.deg).radian, alpha=0.5, zorder=1, color='tab:green', linestyle='--')
       
    if groups is None:
        ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian, sources.dec.wrap_at(180*u.deg).radian,s=size*20)
    else:
        for grp in groups:
            ax.scatter(ax_factor*sources.ra.wrap_at(180*u.deg).radian[grp['Ind']], sources.dec.wrap_at(180*u.deg).radian[grp['Ind']],s=size[grp['Ind']]*20,marker=grp['marker'], color=grp['color'])
    # -------------------------------------------------------------------------
    xt0 = ax_factor*np.array([-135,-90,-45,0,45,90,135])*np.pi/180
    xtl0 =         ['15$^\mathrm{h}$','18$^\mathrm{h}$','21$^\mathrm{h}$','0$^\mathrm{h}$','3$^\mathrm{h}$','6$^\mathrm{h}$','9$^\mathrm{h}$']
    xt1 = ax_factor*np.array([-150,-120,-90,-60,-30,0,30,60,90,120,150])*np.pi/180
    xtl1 =         ['14$^\mathrm{h}$','16$^\mathrm{h}$','18$^\mathrm{h}$','20$^\mathrm{h}$','22$^\mathrm{h}$','0$^\mathrm{h}$','2$^\mathrm{h}$','4$^\mathrm{h}$','6$^\mathrm{h}$','8$^\mathrm{h}$','10$^\mathrm{h}$']
    yt = np.array([-75,-60,-45,-30,-15,0,15,30,45,60,75])*np.pi/180
    ytl = ['-75°','-60°','-45°','-30°','-15°','0°','15°','30°','45°','60°','']
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
    # -------------------------------------------------------------------------
    plt.tight_layout()
    fig.patch.set_facecolor('white')
    if file_out is not None:
        fig.savefig(file_out, dpi=300,format='png',bbox_inches='tight')
    if isshow:
        plt.show()
    plt.close()
    return 0
#===================
def catalog_to_Skycoord(sour_list):
    sour_coor = SkyCoord(ra=sour_list['RA'],dec=sour_list['DEC'],frame='icrs')
    return sour_coor
#===================
def flux_to_groups(flux,steps=[0.5, 1.0, 3.0]):

    groups = []

    colors = ['blue','red','tab:brown','purple']
    markers = ['o','<','s','x']
    grp_lab = ['Weak','Middle','Strong','Insane']

    for i, s in enumerate(steps): # small -> large
        if i == 0:
            ind = np.where(flux < s)
        else:
            ind = np.where((flux >= steps[i-1]) & (flux < steps[i]))
        grp_para = {'Ind':ind,'color':colors[i],'marker':markers[i]}
        groups.append(grp_para)
        print(f'Number of {grp_lab[i]}: {ind[0].shape[0]}')
        if i == len(steps) - 1:
            ind = np.where(flux >= s)
            grp_para = {'Ind':ind,'color':colors[i+1],'marker':markers[i+1]}
            groups.append(grp_para)
            print(f'Number of {grp_lab[i+1]}: {ind[0].shape[0]}')

    return groups
#===================
def main(args):
    # check if file exist
    if not os.path.isfile(args.file_csv):
        print(f"File {args.file_csv} does not exist.")
        return 0

    sour_list = load_catalog_np(args.file_csv)
    sour_coor = catalog_to_Skycoord(sour_list)
    if args.isgroup:
        groups = flux_to_groups(sour_list['Flux'])
    else:
        groups = None
    plot_equ_sky_scatter(sour_coor,size=sour_list['Flux'],fontsize='xx-large', \
            file_out=args.file_out, ax_factor=args.ax_factor, \
            isshow=args.isshow, plotGP=args.plotGP, groups=groups)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
            description='''To plot flux histgram''',
            epilog='version 0.1.')

    parser.add_argument('file_csv', type=str, help='the input catalog csv file') 
    parser.add_argument('--file_out', type=str, help='the output figure name')
    parser.add_argument('--ax_factor', type=int, default=-1,help='X-axis reverse factor, -1')
    parser.add_argument('--plotGP', action='store_true', help='set to add Galactic plane')
    parser.add_argument('--isshow', action='store_true', help='set to show plot')
    parser.add_argument('--isgroup',action='store_true', help='set to plot source in flux groups')

    args = parser.parse_args()
    main(args)

