load('16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N.mat')
load('scales.mat')

% sign is arbitrary, but sign of FPs and tk must be changed together
ndisc = 1; plot_multivariate_model_bias_patterns(-tk,-FPs,ndisc,42,1979,LON_AXIS,LAT_AXIS,sst_scale,slp_scale,linspace(-0.48,0.48,25),linspace(-90,90,25));
ndisc = 2; plot_multivariate_model_bias_patterns(-tk,-FPs,ndisc,42,1979,LON_AXIS,LAT_AXIS,sst_scale,slp_scale,linspace(-0.48,0.48,25),linspace(-90,90,25));