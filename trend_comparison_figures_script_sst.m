
% uncomment line 50 to make individual model figures

label = '1979-2020_Amean'; 
%label = '1958-2021_Amean';

if contains(label,'1958')
    models = {'access','canesm2','canesm5','cesm1','cesm2','csiro_mk36', ...
        'gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
        'miroc_esm2l','mpi','norcpm1'};
    titles = {'ACCESS-ESM1-5','CanESM2','CanESM5','CESM1','CESM2','CSIRO-Mk3-6', ...
        'GFDL-CM3','GFDL-ESM2M','GISS-E2-1-G','IPSL-CM6A-LR','MIROC6', ...
        'MIROC-ES2L','MPI-ESM','NorCPM1'};
else
    models = {'access','canesm2','canesm5','cesm1','cesm2','cnrm_cm6','csiro_mk36', ...
        'ec-earth3','gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
        'miroc_esm2l','mpi','norcpm1'};
    titles = {'ACCESS-ESM1-5','CanESM2','CanESM5','CESM1','CESM2','CNRM-CM6-1','CSIRO-Mk3-6', ...
        'EC-Earth3','GFDL-CM3','GFDL-ESM2M','GISS-E2-1-G','IPSL-CM6A-LR','MIROC6', ...
        'MIROC-ES2L','MPI-ESM','NorCPM1'};
end

varnam = 'tos';  % 'tos', 'psl'
projects_output_dir = '/Users/rcwills/Documents/Data/projects_output_tmp/MMLEA/';
obs = 'ERSST5'; % 'ERSST5', 'COBE', 'AMIPII', 'ERA5', 'JRA55'
load([projects_output_dir,'observations_',varnam,'_',label,'_trend.mat'])
obs_trend = ERSST5_trend; % ERSST5_trend, COBE_trend, AMIPII_trend, ERA5_trend, JRA55_trend

remove_GMSLP = 0; % flag for whether to remove global-mean SLP

ne_all = 0;
trend_std_all = 0;
field1_trends_all = [];
field2_trends_all = [];

for i = 1:length(models)
    if strcmp(models{i},'csiro_mk36') || strcmp(models{i},'gfdl_cm3')
        load([projects_output_dir,models{i},'_ts-sst_psl_',label,'_trends.mat'])
    else
        load([projects_output_dir,models{i},'_tos_psl_',label,'_trends.mat'])
    end
    ne = size(field1_trends,3);
    ne = 10; % uncomment to only use 10-members per model in multi-model analysis (doesn't affect individual model analysis)
    field1_trends_all(:,:,ne_all+1:ne_all+ne) = field1_trends(:,:,1:ne);
    field2_trends_all(:,:,ne_all+1:ne_all+ne) = field2_trends(:,:,1:ne);
    ne_all = ne_all + ne;
    switch varnam
        case 'tos'
            trend_std_all = trend_std_all + ne.*var(field1_trends,0,3);
            %rescaling(i) = trend_comparison_figures(LON_AXIS,LAT_AXIS,field1_trends,obs_trend,linspace(-1.6,1.6,25),linspace(-6,6,25),2,titles{i},'ersst5',label);
        case 'psl'
            %trend_comparison_figures(LON_AXIS,LAT_AXIS,field2_trends,obs_trend,linspace(-240,240,25),linspace(-6,6,25),remove_GMSLP,titles{i},'era5',label);
            if remove_GMSLP==1
                for j = 1:size(field2_trends,3)
                    field2_trends(:,:,j) = field2_trends(:,:,j) - global_mean(LON_AXIS,LAT_AXIS,field2_trends(:,:,j));
                end
            end
            trend_std_all = trend_std_all + ne.*var(field2_trends,0,3);
    end
end

trend_std_all = sqrt(trend_std_all./ne_all);

switch varnam
    case 'tos'
        trend_comparison_figures(LON_AXIS,LAT_AXIS,field1_trends_all,obs_trend,linspace(-1.6,1.6,25),linspace(-5,5,21),2,'Multi-Model Ensemble',obs,label,trend_std_all);
    case 'psl'
        trend_comparison_figures(LON_AXIS,LAT_AXIS,field2_trends_all,obs_trend,linspace(-400,400,25),linspace(-6,6,21),remove_GMSLP,'Multi-Model Ensemble',obs,label,trend_std_all);
end