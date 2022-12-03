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
import numpy as np

def load_sour_list(cata_list,band,sour_hdr):
    with open(cata_list, 'r') as fin:
        i = 0 # line index
        j = 0 # total sour index
        k = 0 # band sour index
        sour_list = []
        for line in fin:
            line_str = line.split()
            if len(line_str) > 1:
                #print(i, j, len(line_str), line_str)
                if line_str[1] == 'J2000':
                    dec_str = line_str[4]
                    dec_str = dec_str.replace("'","m")
                    dec_str = dec_str.replace('"','s')
                    sour = {sour_hdr[0]: k,                # Index
                            sour_hdr[1]: 'J'+line_str[0],  # Name
                            sour_hdr[2]: band,             # Band
                            sour_hdr[3]: line_str[2],      # Position Accuracy
                            sour_hdr[4]: line_str[3],      # RA
                            sour_hdr[5]: dec_str,          # DEC
                            sour_hdr[6]: None,             # Flux 
                            sour_hdr[7]: None,             # Quality
                            }
                    j = j+1
                elif line_str[1] == band and len(line_str) >= 7:
                    #print(i, j, k, len(line_str), line_str)
                    sour['Quality'] = line_str[5] 
                    try:
                        sour['Flux'] = float(line_str[6]) 
                    except ValueError:
                        sour['Flux'] = np.nan 
                    sour_list.append(sour)
                    print(sour)
                    k = k+1
            i = i + 1
    return sour_list

def write_catalog_csv(sour_list, file_out,sour_hdr):
    with open(file_out, 'w') as fout:
        writer = csv.DictWriter(fout, fieldnames=sour_hdr)
        writer.writeheader()
        writer.writerows(sour_list)
    return 0

def main(args):
    # check if file exist
    if not os.path.isfile(args.vla_cata_list):
        print(f"File {args.vla_cata_list} does not exist.")
        return 0

    # check if band is correct:
    if args.band not in ['P', 'L', 'C', 'X', 'U', 'K', 'Q']:
        print(f"Band {args.band} is not correct")
        return 0
    
#    sour_hdr = {'IndexF': -1,'IndexB': -1, 'Name':'J'+line_str[0], 'Band':band, 'PC': line_str[2], 'RA':line_str[3], 'DEC':dec_str, 'Flux':None, 'Quality': None}
    sour_hdr = ['Index', 'Name', 'Band', 'PC', 'RA', 'DEC', 'Flux', 'Quality']
    sour_list = load_sour_list(args.vla_cata_list, args.band,sour_hdr)
    write_catalog_csv(sour_list,args.file_out,sour_hdr)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
            description='''To extract sources and convert to csv from VLA calibrator list''',
            epilog='version 0.1.')

    parser.add_argument('vla_cata_list', type=str, help='the input VLA calibrator list') 
    parser.add_argument('--file_out', type=str, default='cata_out.csv', help='the output catalog file name')
    parser.add_argument('--band', type=str, default='L', help='band to cut, in P, L, C, X, U, K, Q')

    args = parser.parse_args()
    main(args)


