#!/bin/bash -x
mod=${1}

################################################################
# Extract data from the MC packet lists
# Light curves
python src/get_dat.py ${mod}_lc_1-4000_kev 1 4000 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/plt_dat.py lc ${mod}_lc_1-4000_kev ${mod}_lc_1-4000_kev ${mod}_lc_1-4000_kev_min_dir ${mod}_lc_1-4000_kev_max_dir none
python src/plt_dat.py lc ${mod}_lc_1-4000_kev_ele ${mod}_lc_1-4000_kev none

# HEXE
python src/get_dat.py ${mod}_lc_15-45_kev 15 45 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_45-105_kev 45 105 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_105-200_kev 105 200 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_15-200_kev 15 200 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016

python src/plt_dat.py lc ${mod}_lc_15-45_kev ${mod}_lc_15-45_kev ${mod}_lc_15-45_kev_min_dir ${mod}_lc_15-45_kev_max_dir hexe_low
python src/plt_dat.py lc ${mod}_lc_45-105_kev ${mod}_lc_45-105_kev ${mod}_lc_45-105_kev_min_dir ${mod}_lc_45-105_kev_max_dir hexe_mid
python src/plt_dat.py lc ${mod}_lc_45-105_kev_dod_dir ${mod}_lc_45-105_kev ${mod}_lc_45-105_kev_min_dir ${mod}_lc_45-105_kev_max_dir $(find dat -name "${mod}_lc_45-105_kev_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g') none
python src/plt_dat.py lc ${mod}_lc_105-200_kev ${mod}_lc_105-200_kev ${mod}_lc_105-200_kev_min_dir ${mod}_lc_105-200_kev_max_dir hexe_high
python src/plt_dat.py lc ${mod}_lc_15-200_kev ${mod}_lc_15-200_kev ${mod}_lc_15-200_kev_min_dir ${mod}_lc_15-200_kev_max_dir hexe_all
python src/plt_dat.py lc ${mod}_lc_15-45_kev_ele ${mod}_lc_15-45_kev hexe_low
python src/plt_dat.py lc ${mod}_lc_45-105_kev_ele ${mod}_lc_45-105_kev hexe_mid
python src/plt_dat.py lc ${mod}_lc_105-200_kev_ele ${mod}_lc_105-200_kev hexe_high
python src/plt_dat.py lc ${mod}_lc_15-200_kev_ele ${mod}_lc_15-200_kev hexe_all

python src/plt_dat.py lc ${mod}_lc_hexe_all ${mod}_lc_15-45_kev ${mod}_lc_45-105_kev ${mod}_lc_105-200_kev none
python src/plt_dat.py norm ${mod}_lc_hexe_normed ${mod}_lc_15-45_kev ${mod}_lc_45-105_kev ${mod}_lc_105-200_kev none

# Ginga
python src/get_dat.py ${mod}_lc_6-16_kev 6 16 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_16-28_kev 16 28 0 1050 cts lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_6-16_kev_erg 6 16 0 1050 erg lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_lc_16-28_kev_erg 16 28 0 1050 erg lc ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016

python src/plt_dat.py lc ${mod}_lc_6-16_kev ${mod}_lc_6-16_kev_erg ${mod}_lc_6-16_kev_erg_min_dir ${mod}_lc_6-16_kev_erg_max_dir ginga_soft
python src/plt_dat.py lc ${mod}_lc_16-28_kev ${mod}_lc_16-28_kev_erg ${mod}_lc_16-28_kev_erg_min_dir ${mod}_lc_16-28_kev_erg_max_dir ginga_hard

################################################################
# Spectra
python src/get_dat.py ${mod}_spec_185d 1 4e3 170 200 cts spec ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_spec_300d 1 4e3 285 315 cts spec ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_spec_425d 1 4e3 410 440 cts spec ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016

python src/plt_dat.py spec ${mod}_spec_185d ${mod}_spec_185d ${mod}_spec_185d_min_dir ${mod}_spec_185d_max_dir 185
python src/plt_dat.py spec ${mod}_spec_300d ${mod}_spec_300d ${mod}_spec_300d_min_dir ${mod}_spec_300d_max_dir 320
python src/plt_dat.py spec ${mod}_spec_425d ${mod}_spec_425d ${mod}_spec_425d_min_dir ${mod}_spec_425d_max_dir 425

python src/plt_dat.py spec ${mod}_spec_300d_dod_dir ${mod}_spec_300d ${mod}_spec_300d_min_dir ${mod}_spec_300d_max_dir $(find dat -name "${mod}_spec_300d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g') 320

python src/plt_dat.py spec ${mod}_spec_185d_ele ${mod}_spec_185d 185
python src/plt_dat.py spec ${mod}_spec_300d_ele ${mod}_spec_300d 320
python src/plt_dat.py spec ${mod}_spec_425d_ele ${mod}_spec_425d 425

################################################################
#Line light curve
python src/get_dat.py ${mod}_llc847 1 4e3 0 1050 cts llc847 none none ${mod}_015 none
python src/get_dat.py ${mod}_llc1238 1 4e3 0 1050 cts llc1238 none none ${mod}_015 none
python src/get_dat.py ${mod}_llc847+1238 1 4e3 0 1050 cts llc847+1238 none none ${mod}_015 none

python src/plt_dat.py llc847 ${mod}_llc847 ${mod}_llc847 ${mod}_llc847_min_dir ${mod}_llc847_max_dir none
python src/plt_dat.py llc1238 ${mod}_llc1238 ${mod}_llc1238 ${mod}_llc1238_min_dir ${mod}_llc1238_max_dir none
python src/plt_dat.py llc847+1238 ${mod}_llc847+1238 ${mod}_llc847+1238 ${mod}_llc847+1238_min_dir ${mod}_llc847+1238_max_dir none
python src/plt_dat.py llc847+1238 ${mod}_llc847+1238_dod_dir ${mod}_llc847+1238 ${mod}_llc847+1238_min_dir ${mod}_llc847+1238_max_dir $(find dat -name "${mod}_llc847+1238_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g') none

################################################################
# Line profiles
python src/get_dat.py ${mod}_lp847_185d 1 4e3 170 200 cts lp847 none none ${mod}_013 none
python src/get_dat.py ${mod}_lp847_300d 1 4e3 285 315 cts lp847 none none ${mod}_017 none
python src/get_dat.py ${mod}_lp847_425d 1 4e3 410 440 cts lp847 none none ${mod}_014 none
python src/get_dat.py ${mod}_lp1238_185d 1 4e3 170 200 cts lp1238 none none ${mod}_013 none
python src/get_dat.py ${mod}_lp1238_300d 1 4e3 285 315 cts lp1238 none none ${mod}_017 none
python src/get_dat.py ${mod}_lp1238_425d 1 4e3 410 440 cts lp1238 none none ${mod}_014 none

python src/plt_dat.py lp847 ${mod}_lp847 ${mod}_lp847_185d ${mod}_lp847_300d ${mod}_lp847_425d none
python src/plt_dat.py lp847 ${mod}_lp847_185d ${mod}_lp847_185d ${mod}_lp847_185d_min_dir ${mod}_lp847_185d_max_dir none
python src/plt_dat.py lp847 ${mod}_lp847_300d ${mod}_lp847_300d ${mod}_lp847_300d_min_dir ${mod}_lp847_300d_max_dir none
python src/plt_dat.py lp847 ${mod}_lp847_425d ${mod}_lp847_425d ${mod}_lp847_425d_min_dir ${mod}_lp847_425d_max_dir none

python src/plt_dat.py lp847 ${mod}_lp847_300d_dod_dir ${mod}_lp847_300d ${mod}_lp847_300d_min_dir ${mod}_lp847_300d_max_dir $(find dat -name "${mod}_lp847_300d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g') none

python src/plt_dat.py lp1238 ${mod}_lp1238 ${mod}_lp1238_185d ${mod}_lp1238_300d ${mod}_lp1238_425d none
python src/plt_dat.py lp1238 ${mod}_lp1238_185d ${mod}_lp1238_185d ${mod}_lp1238_185d_min_dir ${mod}_lp1238_185d_max_dir none
python src/plt_dat.py lp1238 ${mod}_lp1238_300d ${mod}_lp1238_300d ${mod}_lp1238_300d_min_dir ${mod}_lp1238_300d_max_dir none
python src/plt_dat.py lp1238 ${mod}_lp1238_425d ${mod}_lp1238_425d ${mod}_lp1238_425d_min_dir ${mod}_lp1238_425d_max_dir none

python src/plt_dat.py spec ${mod}_spec_all ${mod}_spec_185d ${mod}_spec_300d ${mod}_spec_425d none
python src/plt_dat.py lp1238 ${mod}_lp ${mod}_lp847_185d ${mod}_lp847_300d ${mod}_lp847_425d ${mod}_lp1238_300d ${mod}_lp847_300d_min_dir ${mod}_lp847_300d_max_dir disable

################################################################
# Maps
python src/get_dat.py ${mod}_map_185d 1 4000 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_185d_15-45_kev 15 45 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_185d_45-105_kev 45 105 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_185d_837-857_kev 832 862 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_185d_1228-1248_kev 1218 1258 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_185d_400-4000_kev 400 4000 170 200 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016

python src/get_dat.py ${mod}_map_300d 1 4000 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_300d_15-45_kev 15 45 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_300d_45-105_kev 45 105 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_300d_837-857_kev 832 862 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_300d_1228-1248_kev 1218 1258 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_300d_400-4000_kev 400 4000 285 315 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016

python src/get_dat.py ${mod}_map_425d 1 4000 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_425d_15-45_kev 15 45 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_425d_45-105_kev 45 105 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_425d_837-857_kev 832 862 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_425d_1228-1248_kev 1218 1258 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_425d_400-4000_kev 400 4000 410 440 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016

python src/plt_map.py ${mod}_map_185d ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_185d_400-4000_kev_smo ${mod} feconix min max

python src/plt_map.py ${mod}_map_300d ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_300d_400-4000_kev_smo ${mod} feconix min max

python src/plt_map.py ${mod}_map_425d ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_425d_400-4000_kev_smo ${mod} feconix min max
