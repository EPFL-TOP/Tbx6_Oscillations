# Tbx6 Oscillations

The repository contains:
- The folder `MastodonTracks` with CSV files generated by the Fiji plugin Mastodon and containing the data from cell tracking used in Fig. 1, 2 and Extended Data Fig. 2, 3.
- The folder `ROI_profiles` with ROI profiles from Maximum-intensity projections of elliptically transformed Tbx6-mNG timelapses and combined immunostainings against Tbx6 and in situ hybdridization for ripply1 mRNA used in Fig. 3, Extended Data Fig. 4. 
- The folder `Python_scripts` with the scripts `Venzinetal_figures.py` that generates all figures and va_utils.py that contains some utility functions. 

The script `Venzinetal_figures.py` should be run with Python > 3.7 and the packages version specified in the `env.yml`

Create the environment and activate it
```shell
conda env create -f env.yml
conda activate Tbx6_Oscillations
```

then go in the folder `Data_Venzin_et_al/Python_scripts` and run 

```shell
python Venzinetal_figures.py
```

The figures will appear in the previous folder `Data_Venzin_et_al`.