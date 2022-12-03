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
import csv

def load_catalog_csv(file_cata):
    sour_list = []
    with open(file_cata,'r') as filein:
        reader = csv.DictReader(filein)
        hdr = reader.fieldnames
        for row in reader:
            sour = {hdr[0]: int(row[hdr[0]])  , # Index
                    hdr[1]: row[hdr[1]]       , # Name
                    hdr[2]: row[hdr[2]]       , # Band
                    hdr[3]: row[hdr[3]]       , # Position Accuracy
                    hdr[4]: row[hdr[4]]       , # RA
                    hdr[5]: row[hdr[5]]       , # DEC
                    hdr[6]: float(row[hdr[6]]), # Flux 
                    hdr[7]: row[hdr[7]]       , # Quality
                    }
            sour_list.append(sour)
    return sour_list, hdr

def write_catalog_csv(sour_list, file_out, sour_hdr):
    with open(file_out, 'w') as fout:
        writer = csv.DictWriter(fout, fieldnames=sour_hdr)
        writer.writeheader()
        writer.writerows(sour_list)
    return 0

def catalog_cut_flux(sour_list, flux):
    sour_list_cut = []
    for sour in sour_list:
        if sour['Flux'] >= flux:
            sour_list_cut.append(sour)
    return sour_list_cut

def main(args):
    # check if file exist
    if not os.path.isfile(args.catalog):
        print(f"File {args.catalog} does not exist.")
        return 0

    # check if band is correct:
    if args.flux < 0:
        print(f"Flux {args.band} must be > 0.")
        return 0
    
    sour_list, sour_hdr = load_catalog_csv(args.catalog)
    sour_list_cut = catalog_cut_flux(sour_list, args.flux)
    write_catalog_csv(sour_list_cut,args.file_out,sour_hdr)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
            description='''To extract sources and convert to csv from VLA calibrator list''',
            epilog='version 0.1.')

    parser.add_argument('catalog', type=str, help='the input VLA calibrator list') 
    parser.add_argument('--file_out', type=str, default='cata_out.csv', help='the output catalog file name')
    parser.add_argument('--flux', type=float, default=1, help='flux threshold to cut, default 1 Jy')
    args = parser.parse_args()
    main(args)


