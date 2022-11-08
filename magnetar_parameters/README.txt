################################################################
run with "h5fc"
h5fc: aliased to /usr/local/Cellar/hdf5/1.10.2_1/bin/h5fc

this is the only way I found to make the compiler find the required library files for hdf5


################################################################
Folders:
'src' contains the final version of the Fortran code and its visualizer

'out' contains output from the Fortran code, i.e. the mc packets. each
data set (r??.esc) should be identified by a two-digit
number and should have an associated meta file (the .txt file)

'dat' contains the intermediate step, the simulation output has been
prepared but is not combined into final products. This means light
curves, spectra, and maps in text, i.e. data begind the figures. These
also includes the information related to which magnetar parameters
were used given by the identifier m??

'fig' contains the final graphical products


################################################################
Workflow:
1. Create SN (toy) model with toy.py
2. Use grt.sh to prepare grt.f and run it, this is the Fortran step
that produces the MC packets.
3. Use run_all.sh to prepare data products from the MC packets. The
first step is to take the packets and compute light curves, spectra,
and maps (dont by get_dat.py). Then make products using plt_dat.py and
plt_map.py


################################################################
Labels:
Models are given by model name followed by s??, where ?? is an index
to allow for several versions of the same model.

Output MC packet files have an additional tag r??, which denotes
Fortran run number. This allows for several MC runs of the same model
to provide better resolution in certain parts of the parameter space.

The m?? tag separates different magnetar parameters.

