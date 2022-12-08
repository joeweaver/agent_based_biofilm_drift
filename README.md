# Quantifying Drift-Fitness Balance Using an Agent-Based Biofilm Model of Identical Heterotrophs Under Low Nutrient Conditions

This repository contains the analysis code supporting the journal article (in submission) "Quantifying Drift-Fitness Balance Using an Agent-Based Biofilm Model of Identical Heterotrophs Under Low Nutrient Conditions"

# Data analysis

The relevant analyses are encapsulated in one R project within the `analyses` folder. The scripts are numbered in order of intended exection. Hard dependencies between scripts are included with commented preconditions.
Two important things for first-time runs:
* `0_download_sim_results.R` need only be run once. It retrieves the simulation results from the project [OSF data repository](https://osf.io/fch3z/).
* `6_GAM_and_comparison.R` uses models which take some time to fit. The results are cached as RDS files. To generate these, some code needs to be uncommented for this first run, and this is documented within the script.

Results from the analysis are generated under the `output` folder which is created as the scripts are run. There are subfolders sorting some results by type. Notably, `si` and `eda` respectiviely contain data presented in Supporting Information and used during exploratory data analysis. The latter is probably not interesing to most but contains details on things such as individual sigmoid fits. Images are saved as high resolution (330 dpi) pngs, tiffs, and pdfs.

All runs are logged in the `log` directory and their names include a timestamp. In general, the logs contain the execution environment info and record which files were generated, along with MD5 checksums.

# Simulations
While we expect most people will just download the results, we haved provide a link to the exact variant of NUFEB used  and the Snakemake-based workflows used to run the simulations.

## Getting and running NUFEB
The NUFEB variant used here will be deposited in the `dev-compute-vol` branch of the [NUFEB-dev](https://github.com/nufeb/NUFEB-dev) repository. *Important* this variant is well-out of date. **We recommend anyone interested in developing new NUFEB simulations use the main branch.** This branch exists for exact reproducibility. The initial simulations were started over a year ago and, to maintain comparability between runs,  the code was not updated. As such, many of the input files are slightly different than modern NUFEB files. Additionally, the build process is, unfortunately, still very complicated.

To build NUFEB, you must use the make commands enabling MPI and USER-NUFEB. The simulations run on the cluster do not require VTK. Please refer to the NUFEB repository for detailed build instructions. We recognize that building NUFEB, especially older versions, is non-trivial and are happy to help.

## Managing the simulation workflows
Two Snakemake-based workflows exist in `simulations/snakemakes`  For any crowding condition, such as 3x3 initial organisms spaced 5 diameters apart, you must first generate a small set of simulations (the 'seed' runs)  with one workflow, then run a second workflow for the parameter sweep.
We have included template files for both types of runs. 

For the example of 3x3 and 5 diameters:
First:
```
cp -r seeds_template_mxn_y_seeds/ 5x5_2.5_seeds
cd 5x5_2.5_seeds
snakemake --profile slurm --latency-wait 60 --rerun-incomplete
```

The exact snakemake parameters may vary based on your computing environment.

Then, after those runs complete:
```
cp -r kinetic_template_mxn_y_seeds/ 5x5_2.5_kinetic
cd 5x5_2.5_kinetic
snakemake --profile slurm --latency-wait 60 --rerun-incomplete
```

We highly suggest running within a `tmux` environment or use a similar mechanism which allows snakemake to persist after logging out.

The most important result is `sweep_colony_outcomes.csv` which is generated in the `results/<NxN_spacing>_default_mu_ks_yield_conc` folder for each kinetic sweep. These results summarize the runs and are the files stored on OSF.io and downloaded by the analysis script.

If interested in a 'test run' we also highly suggest reducing the parameter sweep by altering the `sweep_coeffs` list in the kinetic sweep `Snakefile`

# Execution environments

The HPC environment we used incorporated the following modules, some are defaults and may not be needed
```
Currently Loaded Modules:
  1) gompi/2020a                       13) numactl/2.0.13-GCCcore-10.2.0
  2) Szip/2.1.1-GCCcore-9.3.0          14) UCX/1.9.0-GCCcore-10.2.0
  3) HDF5/1.10.6-gompi-2020a           15) libfabric/1.11.0-GCCcore-10.2.0
  4) binutils/2.35-GCCcore-10.2.0      16) OpenMPI/4.0.5-GCC-10.2.0
  5) GCC/10.2.0                        17) bzip2/1.0.8-GCCcore-10.2.0
  6) XZ/5.2.5-GCCcore-10.2.0           18) ncurses/6.2-GCCcore-10.2.0
  7) zlib/1.2.11-GCCcore-10.2.0        19) libreadline/8.0-GCCcore-10.2.0
  8) libxml2/2.9.10-GCCcore-10.2.0     20) Tcl/8.6.10-GCCcore-10.2.0
  9) libpciaccess/0.16-GCCcore-10.2.0  21) SQLite/3.33.0-GCCcore-10.2.0
 10) hwloc/2.2.0-GCCcore-10.2.0        22) GMP/6.2.0-GCCcore-10.2.0
 11) libevent/2.1.12-GCCcore-10.2.0    23) libffi/3.3-GCCcore-10.2.0
 12) GCCcore/10.2.0                    24) Python/3.8.6-GCCcore-10.2.0
```

Similarly, the R environment used included the following:
```
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] tools     stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stringr_1.4.0      viridis_0.6.2      viridisLite_0.4.0  patchwork_1.1.1    cowplot_1.1.1      tidymv_3.3.2      
 [7] itsadug_2.4.1      plotfunctions_1.4  mgcv_1.8-40        nlme_3.1-158       purrr_0.3.4        ggtext_0.1.1      
[13] ggpubr_0.4.0       ggpmisc_0.5.1      ggpp_0.4.5         broom_1.0.0        latex2exp_0.9.4    ggh4x_0.2.2       
[19] magrittr_2.0.3     logger_0.2.2       gt_0.6.0           tidyr_1.2.0        RColorBrewer_1.1-3 readr_2.1.2       
[25] here_1.0.1         dplyr_1.0.9        osfr_0.2.8         gifski_1.6.6-1     transformr_0.1.4   gganimate_1.0.8   
[31] ggplot2_3.3.6     

loaded via a namespace (and not attached):
 [1] fs_1.5.2           sf_1.0-9           progress_1.2.2     httr_1.4.3         rprojroot_2.0.3    backports_1.4.1   
 [7] utf8_1.2.2         R6_2.5.1           KernSmooth_2.23-20 DBI_1.1.3          colorspace_2.0-3   withr_2.5.0       
[13] gridExtra_2.3      tidyselect_1.1.2   prettyunits_1.1.1  curl_4.3.2         compiler_4.2.0     cli_3.4.1         
[19] quantreg_5.94      SparseM_1.81       xml2_1.3.3         scales_1.2.1       classInt_0.4-8     proxy_0.4-27      
[25] digest_0.6.29      pkgconfig_2.0.3    htmltools_0.5.2    fastmap_1.1.0      rlang_1.0.4        rstudioapi_0.14   
[31] httpcode_0.3.0     farver_2.1.1       generics_0.1.3     jsonlite_1.8.0     car_3.1-0          Matrix_1.5-1      
[37] Rcpp_1.0.9         munsell_0.5.0      fansi_1.0.3        abind_1.4-5        lifecycle_1.0.1    stringi_1.7.8     
[43] carData_3.0-5      MASS_7.3-57        grid_4.2.0         crayon_1.5.1       lattice_0.20-45    splines_4.2.0     
[49] gridtext_0.1.4     hms_1.1.1          pillar_1.7.0       ggsignif_0.6.3     lpSolve_5.6.17     crul_1.3          
[55] glue_1.6.2         vctrs_0.4.1        tzdb_0.3.0         tweenr_2.0.2       MatrixModels_0.5-1 gtable_0.3.0      
[61] assertthat_0.2.1   cachem_1.0.6       e1071_1.7-12       rstatix_0.7.0      class_7.3-20       survival_3.3-1    
[67] tibble_3.1.7       memoise_2.0.1      units_0.8-0        ellipsis_0.3.2   
```



