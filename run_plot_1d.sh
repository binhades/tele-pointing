#!/bin/bash
Help()
{
    echo "Set '-n id ' to choose plot taskes:"
    echo "      'n = 111' To plot catalog histgram"
    echo "      'n = 112' for lam=1e5"
    echo "      'n = 113' for lam=1e6"
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
#   1.1 To plot catalog histgram
# ---------------------------
#     1.1.1 To plot catalog histgram: version 1: Ku-band 
# ---------------------------
        if [ "$id" == "111" ]; then
            file_cata="./catalog/callist_k.csv"
            file_imag="./pic/dist_flux_k.png"
            ./plot_catalog_hist.py ${file_cata} --file_out ${file_imag}
        fi
# =================================

