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

python src/plt_dat.py lc ${mod}_lc_15-45_kev ${mod}_lc_15-45_kev ${mod}_lc_15-45_kev_min_dir ${mod}_lc_15-45_kev_max_dir hexe_low
python src/plt_dat.py lc ${mod}_lc_45-105_kev ${mod}_lc_45-105_kev ${mod}_lc_45-105_kev_min_dir ${mod}_lc_45-105_kev_max_dir hexe_mid
python src/plt_dat.py lc ${mod}_lc_105-200_kev ${mod}_lc_105-200_kev ${mod}_lc_105-200_kev_min_dir ${mod}_lc_105-200_kev_max_dir hexe_high
python src/plt_dat.py lc ${mod}_lc_15-45_kev_ele ${mod}_lc_15-45_kev hexe_low
python src/plt_dat.py lc ${mod}_lc_45-105_kev_ele ${mod}_lc_45-105_kev hexe_mid
python src/plt_dat.py lc ${mod}_lc_105-200_kev_ele ${mod}_lc_105-200_kev hexe_high

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
python src/get_dat.py ${mod}_spec_20d 1 4e3 17 23 cts spec ${mod}_020 ${mod}_010 ${mod}_011 ${mod}_016
python src/get_dat.py ${mod}_spec_50d 1 4e3 42 58 cts spec ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_spec_100d 1 4e3 90 110 cts spec ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_spec_70d 1 4e3 55 85 cts spec ${mod}_020 ${mod}_021 ${mod}_015 ${mod}_016
python src/get_dat.py ${mod}_spec_125d 1 4e3 115 135 cts spec ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_spec_250d 1 4e3 240 260 cts spec ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016

python src/plt_dat.py spec ${mod}_spec_20d ${mod}_spec_20d ${mod}_spec_20d_min_dir ${mod}_spec_20d_max_dir none
python src/plt_dat.py spec ${mod}_spec_50d ${mod}_spec_50d ${mod}_spec_50d_min_dir ${mod}_spec_50d_max_dir none
python src/plt_dat.py spec ${mod}_spec_70d ${mod}_spec_70d ${mod}_spec_70d_min_dir ${mod}_spec_70d_max_dir none
python src/plt_dat.py spec ${mod}_spec_100d ${mod}_spec_100d ${mod}_spec_100d_min_dir ${mod}_spec_100d_max_dir none
python src/plt_dat.py spec ${mod}_spec_125d ${mod}_spec_125d ${mod}_spec_125d_min_dir ${mod}_spec_125d_max_dir none
python src/plt_dat.py spec ${mod}_spec_250d ${mod}_spec_250d ${mod}_spec_250d_min_dir ${mod}_spec_250d_max_dir none

python src/plt_dat.py spec ${mod}_spec_20d_ele ${mod}_spec_20d none
python src/plt_dat.py spec ${mod}_spec_50d_ele ${mod}_spec_50d none
python src/plt_dat.py spec ${mod}_spec_70d_ele ${mod}_spec_70d none
python src/plt_dat.py spec ${mod}_spec_100d_ele ${mod}_spec_100d none
python src/plt_dat.py spec ${mod}_spec_125d_ele ${mod}_spec_125d none
python src/plt_dat.py spec ${mod}_spec_250d_ele ${mod}_spec_250d none


################################################################
#Line light curve
python src/get_dat.py ${mod}_llc847 1 4e3 0 1050 cts llc847 none none ${mod}_015 none
python src/get_dat.py ${mod}_llc1238 1 4e3 0 1050 cts llc1238 none none ${mod}_015 none
python src/get_dat.py ${mod}_llc847+1238 1 4e3 0 1050 cts llc847+1238 none none ${mod}_015 none

python src/plt_dat.py llc847 ${mod}_llc847 ${mod}_llc847 ${mod}_llc847_min_dir ${mod}_llc847_max_dir none
python src/plt_dat.py llc1238 ${mod}_llc1238 ${mod}_llc1238 ${mod}_llc1238_min_dir ${mod}_llc1238_max_dir none
python src/plt_dat.py llc847+1238 ${mod}_llc847+1238 ${mod}_llc847+1238 ${mod}_llc847+1238_min_dir ${mod}_llc847+1238_max_dir none

################################################################
# Line profiles
python src/get_dat.py ${mod}_lp847_50d 1 4e3 42 58 cts lp847 none none ${mod}_013 none
python src/get_dat.py ${mod}_lp847_125d 1 4e3 115 135 cts lp847 none none ${mod}_017 none
python src/get_dat.py ${mod}_lp847_250d 1 4e3 240 260 cts lp847 none none ${mod}_014 none
python src/get_dat.py ${mod}_lp1238_50d 1 4e3 42 58 cts lp1238 none none ${mod}_013 none
python src/get_dat.py ${mod}_lp1238_125d 1 4e3 115 135 cts lp1238 none none ${mod}_017 none
python src/get_dat.py ${mod}_lp1238_250d 1 4e3 240 260 cts lp1238 none none ${mod}_014 none

python src/plt_dat.py lp847 ${mod}_lp847 ${mod}_lp847_50d ${mod}_lp847_125d ${mod}_lp847_250d none
python src/plt_dat.py lp847 ${mod}_lp847_50d ${mod}_lp847_50d ${mod}_lp847_50d_min_dir ${mod}_lp847_50d_max_dir none
python src/plt_dat.py lp847 ${mod}_lp847_125d ${mod}_lp847_125d ${mod}_lp847_125d_min_dir ${mod}_lp847_125d_max_dir none
python src/plt_dat.py lp847 ${mod}_lp847_250d ${mod}_lp847_250d ${mod}_lp847_250d_min_dir ${mod}_lp847_250d_max_dir none
python src/plt_dat.py lp1238 ${mod}_lp1238 ${mod}_lp1238_50d ${mod}_lp1238_125d ${mod}_lp1238_250d none
python src/plt_dat.py lp1238 ${mod}_lp1238_50d ${mod}_lp1238_50d ${mod}_lp1238_50d_min_dir ${mod}_lp1238_50d_max_dir none
python src/plt_dat.py lp1238 ${mod}_lp1238_125d ${mod}_lp1238_125d ${mod}_lp1238_125d_min_dir ${mod}_lp1238_125d_max_dir none
python src/plt_dat.py lp1238 ${mod}_lp1238_250d ${mod}_lp1238_250d ${mod}_lp1238_250d_min_dir ${mod}_lp1238_250d_max_dir none

python src/plt_dat.py spec ${mod}_spec_all ${mod}_spec_50d ${mod}_spec_125d ${mod}_spec_250d none
python src/plt_dat.py lp1238 ${mod}_lp ${mod}_lp847_50d ${mod}_lp847_125d ${mod}_lp847_250d ${mod}_lp1238_125d ${mod}_lp847_125d_min_dir ${mod}_lp847_125d_max_dir disable

################################################################
# Maps
python src/get_dat.py ${mod}_map_50d 1 4000 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_50d_15-45_kev 15 45 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_50d_45-105_kev 45 105 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_50d_837-857_kev 832 862 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_50d_1228-1248_kev 1218 1258 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016
python src/get_dat.py ${mod}_map_50d_400-4000_kev 400 4000 42 58 cts map ${mod}_020 ${mod}_021 ${mod}_013 ${mod}_016

python src/get_dat.py ${mod}_map_100d 1 4000 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_map_100d_15-45_kev 15 45 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_map_100d_45-105_kev 45 105 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_map_100d_837-857_kev 832 862 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_map_100d_1228-1248_kev 1218 1258 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016
python src/get_dat.py ${mod}_map_100d_400-4000_kev 400 4000 90 110 cts map ${mod}_020 ${mod}_021 ${mod}_018 ${mod}_016

python src/get_dat.py ${mod}_map_125d 1 4000 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_125d_15-45_kev 15 45 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_125d_45-105_kev 45 105 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_125d_837-857_kev 832 862 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_125d_1228-1248_kev 1218 1258 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016
python src/get_dat.py ${mod}_map_125d_400-4000_kev 400 4000 115 135 cts map ${mod}_020 ${mod}_021 ${mod}_017 ${mod}_016

python src/get_dat.py ${mod}_map_250d 1 4000 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_250d_15-45_kev 15 45 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_250d_45-105_kev 45 105 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_250d_837-857_kev 832 862 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_250d_1228-1248_kev 1218 1258 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016
python src/get_dat.py ${mod}_map_250d_400-4000_kev 400 4000 240 260 cts map ${mod}_020 ${mod}_021 ${mod}_014 ${mod}_016

python src/plt_map.py ${mod}_map_50d ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_50d_400-4000_kev_smo ${mod} feconix min max

python src/plt_map.py ${mod}_map_100d ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_100d_400-4000_kev_smo ${mod} feconix min max

python src/plt_map.py ${mod}_map_125d ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_125d_400-4000_kev_smo ${mod} feconix min max

python src/plt_map.py ${mod}_map_250d ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_15-45_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_15-45_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_45-105_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_45-105_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_837-857_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_837-857_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_1228-1248_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_1228-1248_kev_smo ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_400-4000_kev ${mod} feconix min max
python src/plt_map.py ${mod}_map_250d_400-4000_kev_smo ${mod} feconix min max
