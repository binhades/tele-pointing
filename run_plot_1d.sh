#!/usr/bin/env bash
Help()
{
    echo "Set '-n id ' to choose plot taskes:"
    echo "      'n = 1' To plot catalog related figures"
}

# ---------------------------
# process ARGs

while getopts n:h flag
do
    case "${flag}" in
        n) id=${OPTARG};; 
        h) # display help
            Help
            exit;;
        \?) # incorrect option
            echo "Error: Invalid option, -h for help"
            exit;;
    esac
done

if [ -n "$id" ]; then
    echo "Working on Task: ${id}"
else
    echo "Set -n for Task id"
    Help
    exit
fi
# =================================
# 1. Plot catalog related figures.
# ---------------------------
#   1.1 To plot catalog flux histogram
# ---------------------------
#     1.1.1 To plot catalog flux histogram: version 1: Ku-band 
# ---------------------------
        if [ "$id" == "111" ]; then
            file_cata="./catalog/vla_callist_Ku-band.csv"
            file_imag="./pic/cata_flux_hist_k.png"
            ./plot_catalog_flux_hist.py ${file_cata} --file_out ${file_imag}
        fi
# =================================

