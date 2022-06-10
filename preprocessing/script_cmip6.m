addpath(genpath('/eos15/rcwills/Climate_analysis/base_code/'))
addpath(genpath('/eos15/rcwills/Climate_analysis/projects/MMLEA'))

year1 = 1979; % 1958
year2 = 2020; % 2021

preprocess_cmip6_LE('access','ssp245',[1:10,24,27,29],year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('canesm5','ssp370',1:25,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('ipsl_cm6a','ssp370',[1:10,14],year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc6','ssp585',1:50,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc_esm2l','ssp245',1:30,year1,year2,'Amean','Omon','tos','Amon','psl');

preprocess_cmip6_LE('access','ssp245',[1:10,24,27,29],year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('canesm5','ssp370',1:25,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('ipsl_cm6a','ssp370',[1:10,14],year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc6','ssp585',1:50,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc_esm2l','ssp245',1:30,year1,year2,'DJF','Omon','tos','Amon','psl');

preprocess_cmip6_LE('access','ssp245',[1:10,24,27,29],year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip6_LE('canesm5','ssp370',1:25,year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip6_LE('ipsl_cm6a','ssp370',[1:10,14],year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc6','ssp585',1:50,year1,year2,'JJA','Omon','tos','Amon','psl');
preprocess_cmip6_LE('miroc_esm2l','ssp245',1:30,year1,year2,'JJA','Omon','tos','Amon','psl');

preprocess_cmip6_LE('ec-earth3','ssp585',101:150,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('ec-earth3','ssp585',101:150,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('ec-earth3','ssp585',101:150,year1,year2,'JJA','Omon','tos','Amon','psl');

cesm2_macro1 = [100101 102102 104103 106104 108105 110106 112107 114108 116109 118110];
cesm2_macro2 = cesm2_macro1+1000;
cesm2_members = [123101:123120, 125101:125120, 128101:128120, 130101:130120, cesm2_macro1, cesm2_macro2(1:9)];

preprocess_cmip6_LE('cesm2','ssp370',cesm2_members,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('cesm2','ssp370',cesm2_members,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('cesm2','ssp370',cesm2_members,year1,year2,'JJA','Omon','tos','Amon','psl');

giss_members = 112:100:1012;

preprocess_cmip6_LE('giss_e21g','ssp370',giss_members,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('giss_e21g','ssp370',giss_members,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('giss_e21g','ssp370',giss_members,year1,year2,'JJA','Omon','tos','Amon','psl');

preprocess_cmip6_LE('cnrm_cm6','ssp245',1:10,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('cnrm_cm6','ssp245',1:10,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('cnrm_cm6','ssp245',1:10,year1,year2,'JJA','Omon','tos','Amon','psl');

preprocess_cmip6_LE('norcpm1','historical-ext',1:30,year1,year2,'Amean','Omon','tos','Amon','psl');
preprocess_cmip6_LE('norcpm1','historical-ext',1:30,year1,year2,'DJF','Omon','tos','Amon','psl');
preprocess_cmip6_LE('norcpm1','historical-ext',1:30,year1,year2,'JJA','Omon','tos','Amon','psl');
