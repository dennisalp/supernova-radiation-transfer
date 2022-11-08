#!/bin/bash -x

mag=("2 0.3 2 1" "2 0.3 1 1" "1 0.1 2 1" "1 0.1 1 1")

for ii in "${!mag[@]}"
do
    if [ "$ii" -eq 0 ] || [ "$ii" -eq 2 ] || [ "$ii" -eq 3 ]  # Excludes 98 and 99.
    then
        continue
    fi

    id=$(printf m%02d $ii)
    
    ################################################################
    # Extract data from the MC packet lists
    # Light curves
    
    ################################################################
    # Spectra
    python plt_dat.py spec 000_001_${id}_spec_100d_dod_dir  slsn-i_000_${id}_spec_100d  $(find ../dat/ -name "slsn-i_001_${id}_spec_100d_dir_0??.txt"  -exec basename {} \; | sed 's/\.txt//g')
    python plt_dat.py spec 000_001_${id}_spec_300d_dod_dir  slsn-i_000_${id}_spec_300d  $(find ../dat/ -name "slsn-i_001_${id}_spec_300d_dir_0??.txt"  -exec basename {} \; | sed 's/\.txt//g')
    python plt_dat.py spec 000_001_${id}_spec_3000d_dod_dir slsn-i_000_${id}_spec_3000d $(find ../dat/ -name "slsn-i_001_${id}_spec_3000d_dir_0??.txt" -exec basename {} \; | sed 's/\.txt//g')

    ################################################################
    # Maps
done

################################################################
# Extract data from the MC packet lists
# Light curves

################################################################
# Spectra
python plt_dat.py spec 000_001_mmm_spec_100d_dod_dir  $(find ../dat/ -name "slsn-i_???_m??_spec_100d.txt"   -exec basename {} \; | sed 's/\.txt//g')
python plt_dat.py spec 000_001_mmm_spec_300d_dod_dir  $(find ../dat/ -name "slsn-i_???_m??_spec_300d.txt"   -exec basename {} \; | sed 's/\.txt//g')
python plt_dat.py spec 000_001_mmm_spec_3000d_dod_dir $(find ../dat/ -name "slsn-i_???_m??_spec_3000d.txt"  -exec basename {} \; | sed 's/\.txt//g')

################################################################
# Maps
