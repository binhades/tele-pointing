#!/bin/bash
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
#   1.1 To plot catalog sky distribution
# ---------------------------
#     1.1.1 To plot: Ku-band Full catalog
# ---------------------------
        if [ "$id" == "111" ]; then
            file_cata="./catalog/vla_callist_Ku-band.csv"
            file_imag="./pic/cata_sky_dist_Ku.png"
            ./plot_catalog_sky_dist.py ${file_cata} --file_out ${file_imag} \
                --isgroup --plotGP --isshow
        fi
# ---------------------------
#     1.1.2 To plot: Ku-band Flux >= 3 catalog
# ---------------------------
        if [ "$id" == "112" ]; then
            file_cata="./catalog/vla_callist_Ku-band_3Jy.csv"
            file_imag="./pic/cata_sky_dist_Ku_3Jy.png"
            ./plot_catalog_sky_dist.py ${file_cata} --file_out ${file_imag} --plotGP --isshow
        fi
# ---------------------------
#     1.1.2 To plot: Ku-band Flux >= 3 catalog
# ---------------------------
        if [ "$id" == "113" ]; then
            file_cata="./catalog/vla_callist_Ku-band_1Jy.csv"
            file_imag="./pic/cata_sky_dist_Ku_1Jy.png"
            ./plot_catalog_sky_dist.py ${file_cata} --file_out ${file_imag} --plotGP --isshow
        fi
# =================================

