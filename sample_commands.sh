# Run the Fortran code
h5fc -fopenmp -march=native -O3 src/grt.f; time ./a.out; cat out/tmp*.esc > out/015.esc; cat out/tmp*.abs > out/015.abs

# Extract data from the MC packet lists
python src/get_dat.py b15_lp1238_425d 1 4e3 410 440 cts lp1238 none none 014 none
python src/get_dat.py n20_lp1238_425d 1 4e3 410 440 cts lp1238 none none 114 none
python src/get_dat.py b15_llc1238 1 4e3 0 1000 cts llc1238 none none 015 016

# Visualize the data
python src/plt_dat.py lc b15_lc_15-45_kev b15_lc_15-45_kev_min_dir b15_lc_15-45_kev_max_dir hexe_low
python src/plt_dat.py spec b15_spec_425d b15_spec_425d_min_dir b15_spec_425d_max_dir 425
python src/plt_map.py b15_map_185d feconix ns min max

python src/plt_dat.py lp1238 b15_lp1238_425d b15_lp1238_425d b15_lp1238_425d_min_dir b15_lp1238_425d_max_dir 425
python src/plt_dat.py lp1238 n20_lp1238_425d n20_lp1238_425d n20_lp1238_425d_min_dir n20_lp1238_425d_max_dir 425
python src/plt_dat.py llc1238 n20_llc1238 n20_llc1238 n20_llc1238_min_dir n20_llc1238_max_dir 425
python src/plt_dat.py llc1238 b15_llc1238 b15_llc1238 b15_llc1238_min_dir b15_llc1238_max_dir 425
