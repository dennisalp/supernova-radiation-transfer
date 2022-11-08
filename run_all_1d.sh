mo1=b15_lmc
mo2=b15_l1d
mo3=hmm

python src/plt_dat.py lc 1D_lc_1-4000_kev ${mo1}_lc_1-4000_kev ${mo2}_lc_1-4000_kev ${mo3}_lc_1-4000_kev none
python src/plt_dat.py lc 1D_lc_15-45_kev ${mo1}_lc_15-45_kev ${mo2}_lc_15-45_kev ${mo3}_lc_15-45_kev hexe_low
python src/plt_dat.py lc 1D_lc_45-105_kev ${mo1}_lc_45-105_kev ${mo2}_lc_45-105_kev ${mo3}_lc_45-105_kev hexe_mid
python src/plt_dat.py lc 1D_lc_105-200_kev ${mo1}_lc_105-200_kev ${mo2}_lc_105-200_kev ${mo3}_lc_105-200_kev hexe_high
python src/plt_dat.py lc 1D_lc_6-16_kev ${mo1}_lc_6-16_kev ${mo2}_lc_6-16_kev ${mo3}_lc_6-16_kev ginga_soft
python src/plt_dat.py lc 1D_lc_16-28_kev ${mo1}_lc_16-28_kev ${mo2}_lc_16-28_kev ${mo3}_lc_16-28_kev ginga_hard
python src/plt_dat.py spec 1D_spec_185d ${mo1}_spec_185d ${mo2}_spec_185d ${mo3}_spec_185d 185
python src/plt_dat.py spec 1D_spec_300d ${mo1}_spec_300d ${mo2}_spec_300d ${mo3}_spec_300d 320
python src/plt_dat.py spec 1D_spec_425d ${mo1}_spec_425d ${mo2}_spec_425d ${mo3}_spec_425d 425
python src/plt_dat.py llc847+1238 1D_llc847+1238 ${mo1}_llc847+1238 ${mo2}_llc847+1238 ${mo3}_llc847+1238 none
python src/plt_dat.py lp847 1D_lp_300d ${mo1}_lp847_300d ${mo2}_lp847_300d ${mo3}_lp847_300d disable
