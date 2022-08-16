%% Parameters

truncation = 20; % number of EOFs 
cutoff = 5; % years
latmax = 80;

varnam = 'multivariate';

%%

label = '1979-2020_Amean';
time = 1979:2020;

if contains(label,'1958')
    models = {'access','canesm2','canesm5','cesm1','cesm2','csiro_mk36', ...
        'gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
        'miroc_esm2l','mpi','norcpm1'};
    
else
    models = {'access','canesm2','canesm5','cesm1','cesm2','cnrm_cm6','csiro_mk36', ...
        'ec-earth3','gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
        'miroc_esm2l','mpi','norcpm1'};
end

projects_output_dir = '/Users/rcwills/Documents/Data/projects_output_tmp/MMLEA/';

if strcmp(varnam,'multivariate')
    obs = 'ersst5';
    load([projects_output_dir,obs,'_',label,'_datamatrix.mat'])
    nt = size(xtf,1);
    field_obs = reshape(insert_cols(xtf,icol_ret,icol_disc)',[144 72 nt]);
    sst_scale = sqrt(global_mean(LON_AXIS,LAT_AXIS,var(field_obs,0,3)));
    field_obs = field_obs./sst_scale;
    if contains(label,'1958')
        obs = 'jra55';
    else
        obs = 'era5';
    end
    load([projects_output_dir,obs,'_',label,'_datamatrix.mat'])
    field = reshape(insert_cols(xtf,icol_ret,icol_disc)',[144 72 nt]);
    slp_scale = sqrt(global_mean(LON_AXIS,LAT_AXIS,var(field,0,3)));
    field_obs(145:288,:,:) = field./slp_scale;
    nlon = 288;
    LON_AXIS = [LON_AXIS; LON_AXIS+360];
    name = ['16model_multivariate_ERSST5_ERA5_',label,'_',num2str(latmax),'deg'];
else
    obs = 'ersst5';
    load([projects_output_dir,obs,'_',label,'_datamatrix.mat'])
    nt = size(xtf,1);
    field_obs = reshape(insert_cols(xtf,icol_ret,icol_disc)',[144 72 nt]);
    nlon = 144;
    name = ['16model_ERSST5_',label,'_',num2str(latmax),'deg'];
end
nlat = 72;

%% Loading

ntotal = 0;
nt10 = 10.*nt;

clear field_all

for i = 1:length(models)
    disp(['Loading model ',num2str(i)])
    switch varnam
        case 'tos'
            if strcmp(models{i},'csiro_mk36') || strcmp(models{i},'gfdl_cm3')
                load([projects_output_dir,models{i},'_ts-sst_psl_',label,'_datamatrix.mat'],'xtf1_all','icol_ret1','icol_disc1')
            else
                load([projects_output_dir,models{i},'_tos_psl_',label,'_datamatrix.mat'],'xtf1_all','icol_ret1','icol_disc1')
            end
            n = size(xtf1_all,1);
            field = reshape(insert_cols(xtf1_all,icol_ret1,icol_disc1)',[144 72 n]);
            field_all(:,:,ntotal+1:ntotal+n) = field;
            is_tenpermodel(ntotal+1:ntotal+nt10) = 1;
            field_tenpermodel(:,:,(i-1)*nt10+1:i*nt10) = reshape(insert_cols(xtf1_all(1:nt10,:),icol_ret1,icol_disc1)',[144 72 nt10]);
            ntotal = ntotal + n;
        case 'psl'
            if strcmp(models{i},'csiro_mk36') || strcmp(models{i},'gfdl_cm3')
                load([projects_output_dir,models{i},'_ts-sst_psl_',label,'_datamatrix.mat'],'xtf2_all','icol_ret2','icol_disc2')
            else
                load([projects_output_dir,models{i},'_tos_psl_',label,'_datamatrix.mat'],'xtf2_all','icol_ret2','icol_disc2')
            end
            n = size(xtf2_all,1);
            field = reshape(insert_cols(xtf2_all,icol_ret2,icol_disc2)',[144 72 n]);
            field_all(:,:,ntotal+1:ntotal+n) = field;
            is_tenpermodel(ntotal+1:ntotal+nt10) = 1;
            field_tenpermodel(:,:,(i-1)*nt10+1:i*nt10) = reshape(insert_cols(xtf2_all(1:nt10,:),icol_ret2,icol_disc2)',[144 72 nt10]);
            ntotal = ntotal + n;
        case 'multivariate'
            if strcmp(models{i},'csiro_mk36') || strcmp(models{i},'gfdl_cm3')
                load([projects_output_dir,models{i},'_ts-sst_psl_',label,'_datamatrix.mat'],'xtf1_all','icol_ret1','icol_disc1')
            else
                load([projects_output_dir,models{i},'_tos_psl_',label,'_datamatrix.mat'],'xtf1_all','icol_ret1','icol_disc1')
            end
            n = size(xtf1_all,1);
            field = reshape(insert_cols(xtf1_all,icol_ret1,icol_disc1)',[144 72 n]);
            field_all(1:144,:,ntotal+1:ntotal+n) = field./sst_scale;
            is_tenpermodel(ntotal+1:ntotal+nt10) = 1;
            field_tenpermodel(1:144,:,(i-1)*nt10+1:i*nt10) = reshape(insert_cols(xtf1_all(1:nt10,:),icol_ret1,icol_disc1)',[144 72 nt10])./sst_scale;
            if strcmp(models{i},'csiro_mk36') || strcmp(models{i},'gfdl_cm3')
                load([projects_output_dir,models{i},'_ts-sst_psl_',label,'_datamatrix.mat'],'xtf2_all','icol_ret2','icol_disc2')
            else
                load([projects_output_dir,models{i},'_tos_psl_',label,'_datamatrix.mat'],'xtf2_all','icol_ret2','icol_disc2')
            end
            field = reshape(insert_cols(xtf2_all,icol_ret2,icol_disc2)',[144 72 n]);
            field_all(145:288,:,ntotal+1:ntotal+n) = field./slp_scale;
            field_tenpermodel(145:288,:,(i-1)*nt10+1:i*nt10) = reshape(insert_cols(xtf2_all(1:nt10,:),icol_ret2,icol_disc2)',[144 72 nt10])./slp_scale;
            ntotal = ntotal + n;
    end 
end

%%

ne = ntotal./nt;

field_all = -(field_all - repmat(field_obs,[1 1 ne]));
field_tenpermodel = -(field_tenpermodel - repmat(field_obs,[1 1 10.*length(models)]));
is_tenpermodel(ntotal) = 0;

field_emean = squeeze(mean(reshape(field_tenpermodel,[nlon nlat nt 10.*length(models)]),4));

%% Preprocessing

[Y,~] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
area(isnan(mean(field_all,3))) = 0; 
area(isnan(mean(field_obs,3))) = 0; 

domain = ones(size(area));

% Can specify sub-domains to analyze using the following form:

domain(Y>latmax) = 0;
domain(Y<-latmax) = 0;

s = size(field_all);
X = reshape(field_all,s(1)*s(2),s(3))'; 
X_tenpermodel = reshape(field_tenpermodel,s(1)*s(2),10.*nt*(length(models)))'; 
Xe = reshape(field_emean,s(1)*s(2),nt)';
%s = size(sst_trends_all);
%X = reshape(sst_trends_all,s(1)*s(2),s(3))';
%Xe = reshape(ERSST_trend,s(1)*s(2),1)';
AREA_WEIGHTS = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret = find(AREA_WEIGHTS~=0 & domain);
icol_disc = find(AREA_WEIGHTS==0 | ~domain);
X = X(:,icol_ret);
X_tenpermodel = X_tenpermodel(:,icol_ret);
Xe = Xe(:,icol_ret);
AREA_WEIGHTS = AREA_WEIGHTS(icol_ret);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS' ./ sum(AREA_WEIGHTS);
scale    = sqrt(normvec);

%% replaced by Xe_f line below
% Xe_f = Xe; t = time';
% for i = 1:size(Xe,2)
%     p = polyfit(t,Xe(:,i),1);
%     tmp = Xe(:,i)-p(1)*t-p(2);
%     tmp = lanczos([flipud(tmp); tmp; flipud(tmp)],1,cutoff);
%     Xe_f(:,i) = tmp((nt+1):2*nt) + (p(1)*t + p(2));
% end

%% Calculate forced patterns / Signal-to-noise maximizing EOF analysis

Xe_f = lanczos_filter_datamatrix_trend_passthrough(Xe,cutoff);

Ct = cov(X_tenpermodel); % only use one ensemble member from each model to compute covariance matrix

[tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs, pcvec, evl, rest] = forced_pattern_analysis(X, Xe_f, truncation, scale,Ct);

save([name,'_S2N_preprocessed.mat'],'X','X_obs','X_tenpermodel','Xe','evl','icol_disc','icol_ret','is_tenpermodel','nt','ntotal','pcvec','rest','scale','time')

% note that 's' is not signal-to-noise ratio, because Ct is computed from
% X_onepermodel while X is used in forced_pattern_analysis

FPs       = insert_cols(FPs, icol_ret, icol_disc);
EOFs       = insert_cols(EOFs, icol_ret, icol_disc);
fingerprintsf        = insert_rows(fingerprints, icol_ret, icol_disc);
