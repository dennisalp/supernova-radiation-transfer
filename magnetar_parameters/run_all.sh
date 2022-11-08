#!/bin/bash -x
mod=${1}
run=${2}

mag=("2 0.3 2 1" "2 0.3 1 1" "1 0.1 2 1" "1 0.1 1 1")

for ii in "${!mag[@]}"
do
    if [ "$ii" -eq 99 ] || [ "$ii" -eq 99 ] || [ "$ii" -eq 98 ]  # Excludes 98 and 99.
    then
        continue
    fi
    
    id=$(printf m%02d $ii)
    dat=${mod}_${run}
    lab=${mod}_${id}
    
    ################################################################
    # Extract data from the MC packet lists
    # Light curves
    # python get_dat.py ${lab}_lc_1-10000_kev 1 1e4 0 10000 lc ${dat} ${mag[$ii]} erg
    # python get_dat.py ${lab}_lc_3-78_kev    3  78 0 10000 lc ${dat} ${mag[$ii]} erg
    
    # python plt_dat.py lc ${lab}_lc ${lab}_lc_1-10000_kev ${lab}_lc_3-78_kev
    # python plt_dat.py lc ${lab}_lc_1-10000_kev_dod_dir ${lab}_lc_1-10000_kev $(find ../dat/ -name "${lab}_lc_1-10000_kev_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')
    # python plt_dat.py lc ${lab}_lc_3-78_kev_dod_dir    ${lab}_lc_3-78_kev    $(find ../dat/ -name "${lab}_lc_3-78_kev_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')
    
    # ################################################################
    # # Spectra
    # python get_dat.py ${lab}_spec_50d    1 1e4   45    55 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_70d    1 1e4   63    77 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_100d   1 1e4   90   110 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_150d   1 1e4  135   155 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_300d   1 1e4  270   330 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_1000d  1 1e4  970  1030 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_3000d  1 1e4 2700  3300 spec ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_spec_10000d 1 1e4 9800 10200 spec ${dat} ${mag[$ii]} cts
    
    # python plt_dat.py spec ${lab}_spec ${lab}_spec_50d ${lab}_spec_70d ${lab}_spec_100d ${lab}_spec_150d ${lab}_spec_300d ${lab}_spec_1000d ${lab}_spec_3000d ${lab}_spec_10000d
    # python plt_dat.py spec ${lab}_spec_100d_dod_dir ${lab}_spec_100d $(find ../dat/ -name "${lab}_spec_100d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')
    # python plt_dat.py spec ${lab}_spec_300d_dod_dir ${lab}_spec_300d $(find ../dat/ -name "${lab}_spec_300d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')
    # python plt_dat.py spec ${lab}_spec_3000d_dod_dir ${lab}_spec_3000d $(find ../dat/ -name "${lab}_spec_3000d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')
    
    # ################################################################
    # # Maps
    # python get_dat.py ${lab}_map_50d           1 4000 42 58 map ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_map_50d_3-78_kev  3   78 42 58 map ${dat} ${mag[$ii]} cts
    
    # python get_dat.py ${lab}_map_125d          1 4000 115 135 map ${dat} ${mag[$ii]} cts
    # python get_dat.py ${lab}_map_125d_3-78_kev 3   78 115 135 map ${dat} ${mag[$ii]} cts
    
    # python plt_map.py ${lab}_map_50d ${lab}
    # python plt_map.py ${lab}_map_50d_smo ${lab}
    # python plt_map.py ${lab}_map_50d_3-78_kev ${lab}
    # python plt_map.py ${lab}_map_50d_3-78_kev_smo ${lab}
    
    # python plt_map.py ${lab}_map_125d ${lab}
    # python plt_map.py ${lab}_map_125d_smo ${lab}
    # python plt_map.py ${lab}_map_125d_3-78_kev ${lab}
    # python plt_map.py ${lab}_map_125d_3-78_kev_smo ${lab}

    python get_dat.py ${lab}_map_3000d          1 4000 2700 3300 map ${dat} ${mag[$ii]} cts
    python get_dat.py ${lab}_map_3000d_3-78_kev 3   78 2700 3300 map ${dat} ${mag[$ii]} cts
    python plt_map.py ${lab}_map_3000d ${lab}
    python plt_map.py ${lab}_map_3000d_smo ${lab}
    python plt_map.py ${lab}_map_3000d_3-78_kev ${lab}
    python plt_map.py ${lab}_map_3000d_3-78_kev_smo ${lab}
done
