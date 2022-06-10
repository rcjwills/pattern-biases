addpath(genpath('/eos15/rcwills/Climate_analysis/base_code/'))
addpath(genpath('/eos15/rcwills/Climate_analysis/projects/MMLEA'))

year1 = 1979; % 1958
year2 = 2020; % 2021

preprocess_cmip5_LE('canesm2',year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip5_LE('cesm1',year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip5_LE('gfdl_esm2m',year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip5_LE('mpi',year1,year2,'Amean','Omon','tos','Amon','psl');

preprocess_cmip5_LE('canesm2',year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip5_LE('cesm1',year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip5_LE('gfdl_esm2m',year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip5_LE('mpi',year1,year2,'DJF','Omon','tos','Amon','psl');

preprocess_cmip5_LE('canesm2',year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip5_LE('cesm1',year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip5_LE('gfdl_esm2m',year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip5_LE('mpi',year1,year2,'JJA','Omon','tos','Amon','psl');

preprocess_cmip5_LE('gfdl_cm3',year1,year2,'Amean','Amon','ts-sst','Amon','psl');
preprocess_cmip5_LE('gfdl_cm3',year1,year2,'DJF','Amon','ts-sst','Amon','psl');
preprocess_cmip5_LE('gfdl_cm3',year1,year2,'JJA','Amon','ts-sst','Amon','psl');

using ts-sst for CSIRO because tos is NaN under sea ice
preprocess_cmip5_LE('csiro_mk36',year1,year2,'Amean','Amon','ts-sst','Amon','psl');
preprocess_cmip5_LE('csiro_mk36',year1,year2,'DJF','Amon','ts-sst','Amon','psl');
preprocess_cmip5_LE('csiro_mk36',year1,year2,'JJA','Amon','ts-sst','Amon','psl');