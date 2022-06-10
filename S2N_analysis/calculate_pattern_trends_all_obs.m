
% create multi-model mean datamatrix
load('16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_S2N_preprocessed.mat')
Xmmm = reshape(X_tenpermodel',[14870 42 160]);
Xmmm = -(Xmmm - repmat(X_obs',[1 1 160]));
Xmmm = mean(Xmmm,3);
Xmmm = insert_cols(Xmmm',icol_ret,icol_disc);

sst_names = {'ersst5','cobe','amip2'};
slp_names = {'era5','jra55'};

load('scales.mat')
for i = 1:length(sst_names)
    load([sst_names{i},'_1979-2020_Amean_datamatrix.mat'])
    X1 = reshape(insert_cols(xtf,icol_ret,icol_disc)',[144 72 42]);
    for j = 1:length(slp_names)
        load([slp_names{j},'_1979-2020_Amean_datamatrix.mat'])
        X2 = reshape(xtf',[144 72 42]);
        
        % create obs datamatrix
        X_newobs = [X1./sst_scale; X2./slp_scale];
        X_newobs = reshape(X_newobs,[288*72 42])';
        
        xtf = X_newobs-Xmmm;
        xtf(isnan(xtf)) = 0;
        
        load('16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N.mat', 'fingerprintsf')
        fingerprintsf(isnan(fingerprintsf)) = 0;
        tk_obs(:,:,i,j) = xtf*fingerprintsf;
    end
end

save('/Users/rcwills/Documents/Data/projects_output_tmp/MMLEA/S2N_Patterns/tk_obs.mat','tk_obs','sst_names','slp_names')