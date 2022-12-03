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

def load_catalog_np(file_cata):
    sour_list = np.genfromtxt(file_cata,delimiter=',',dtype=None,names=True,encoding=None)
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

def main(args):
    # check if file exist
    if not os.path.isfile(args.file_csv):
        print(f"File {args.file_csv} does not exist.")
        return 0

    sour_list = load_catalog_np(args.file_csv)
    plot_hist(sour_list['Flux'],file_out=args.file_out)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
            description='''To plot flux histgram''',
            epilog='version 0.1.')

    parser.add_argument('file_csv', type=str, help='the input catalog csv file') 
    parser.add_argument('--file_out', type=str, help='the output figure name')

    args = parser.parse_args()
    main(args)

