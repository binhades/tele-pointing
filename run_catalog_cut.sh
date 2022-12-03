#!/bin/bash
Help()
{
    echo "Set '-n id ' to choose plot taskes:"
    echo "      'n = 1' To input raw catalog"
    echo "      'n = 2' To input Ku-band full catalog"
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
# 1. Input from the raw vla catalog "txt".
# ---------------------------
#   1.1 Cut for Frequency band
# ---------------------------
#     1.1.1 Ku-band 
# ---------------------------
        if [ "$id" == "111" ]; then
            file_cata_inp="./catalog/vla_callist.txt"
            file_cata_out="./catalog/vla_callist_Ku-band.csv"
            band='U'
            ./catalog_cut-band_raw2csv.py ${file_cata_inp} --file_out ${file_cata_out} \
                --band ${band}
        fi
# =================================
# 2. Input from the Ku-band full catalog "csv".
# ---------------------------
#   2.1 Cut for Flux above source
# ---------------------------
#     2.1.1 Flux > 3 
# ---------------------------
        if [ "$id" == "211" ]; then
            flux=3
            file_cata_inp="./catalog/vla_callist_Ku-band.csv"
            file_cata_out="./catalog/vla_callist_Ku-band_${flux}Jy.csv"
            ./catalog_cut-flux_csv.py ${file_cata_inp} --file_out ${file_cata_out} \
                --flux ${flux}
        fi
# ---------------------------
#     2.1.2 Flux > 1 
# ---------------------------
        if [ "$id" == "212" ]; then
            flux=1
            file_cata_inp="./catalog/vla_callist_Ku-band.csv"
            file_cata_out="./catalog/vla_callist_Ku-band_${flux}Jy.csv"
            ./catalog_cut-flux_csv.py ${file_cata_inp} --file_out ${file_cata_out} \
                --flux ${flux}
        fi
# =================================

