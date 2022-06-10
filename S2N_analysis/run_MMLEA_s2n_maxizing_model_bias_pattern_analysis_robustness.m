%% Parameters

truncation = 20; % number of EOFs 
cutoff = 5; % years
latmax = 80;

savedir = '/eos15/rcwills/projects_output/MMLEA/S2N_Patterns/';
%savedir = '/Users/rcwills/Documents/Data/projects_output_tmp/MMLEA/S2N_Patterns/';
name = '16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg';

load([name,'_S2N_preprocessed.mat'])

models = {'access','canesm2','canesm5','cesm1','cesm2','cnrm_cm6','csiro_mk36', ...
    'ec-earth3','gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
    'miroc_esm2l','mpi','norcpm1'};

X_tenpermodel = X_tenpermodel - repmat(X_obs,[10.*length(models) 1]); %remove field_obs from field_tenpermodel, to be replaced with individual ensemble members
X = X - repmat(X_obs,[598 1]);
% note original sign of both was field_obs - field_tenpermodel and field_obs - field_all

%% Robustness from random sampling within ensembles

ne_all = [13 50 25 40 99 10 30 50 20 30 10 11 50 30 100 30];
cum_ne = cumsum(ne_all);
nmodels = length(ne_all);

if(0)

n_resample = 30; % for bootstrapping, total
n_sample = 10*length(models); % size of each bootstrapped sample

% for comparing with single-model analysis:
% n_sample = 30; % to test dependence on sample size

for i = 1:n_resample
    Xi = [];
    for j = 1:n_sample % to make Xi same size as X_tenpermodel
        model1 = randi(length(models));
        if model1>1
            member1 = cum_ne(model1-1)+randi(ne_all(model1));
        else
            member1 = randi(ne_all(model1));
        end
        is1 = (member1-1)*nt+1:member1*nt;
        member2 = member1;
        while member2 == member1 % to avoid choosing same member
            model2 = randi(length(models));
            if model2>1
                member2 = cum_ne(model2-1)+randi(ne_all(model2));
            else
                member2 = randi(ne_all(model2));
            end
        end
        is2 = (member2-1)*nt+1:member2*nt;
        Xi = [Xi; X(is1,:)-X(is2,:)];
    end
    
    Xi_emean = squeeze(mean(reshape(Xi,[nt n_sample size(Xi,2)]),2));
    Xi_emean_f = lanczos_filter_datamatrix_trend_passthrough(Xi_emean,cutoff);
    
    Ct = cov(Xi);
    [tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs, pcvec, evl, rest] = forced_pattern_analysis(Xi, Xi_emean_f, truncation, scale,Ct);
    FPs       = insert_cols(FPs, icol_ret, icol_disc);
    tk_random(:,:,i) = tk(:,1:2);
    FPs_random(:,:,i) = FPs(1:2,:);
    s_random(:,i) = s;
end

save([savedir,name,'_',num2str(cutoff),'yr_N',num2str(truncation),'_S2N_robustness_random_crossmodel_sampling.mat'],'tk_random','FPs_random','s_random','n_sample')

end

%% Robustness from sampling within ensembles

n_resample = 5; % per model

istart = 0;
nt10 = 10.*nt;

for i = 1:nmodels
    clear tk_single tk_other FPs_single FPs_other s_single s_other
    nei = ne_all(i);
    is = istart+1:istart+nei*nt;
    switch i
        case 1
            is_other = istart+nei*nt+1:ntotal;
            is_other_ten = i*nt10+1:length(models)*nt10;
        case 16
            is_other = 1:istart;
            is_other_ten = 1:(i-1)*nt10;
        otherwise
            is_other = [1:istart, istart+nei*nt+1:ntotal];
            is_other_ten = [1:(i-1)*nt10, i*nt10+1:length(models)*nt10];
    end
    istart = istart+nei*nt;
    X_tenpermodel_other = X_tenpermodel(is_other_ten,:);
    is_tenpermodel_other = is_tenpermodel(is_other_ten);
    X_other = X(is_other,:);
    X_single_model = X(is,:);
    randomize_members = randperm(nei); members = randomize_members(1:n_resample);
    nother = size(X_other,1)./nt;
    
    for j = 1:length(members)
        js = (members(j)-1)*nt+1:members(j)*nt;
        X_single_member = X_single_model(js,:);
        js_other = [1:js(1)-1, js(end)+1:nt*nei];
        X_other_members = X_single_model(js_other,:);
                
        % construct data matrices to input to analysis (based on other
        % models)
        Xi = repmat(X_single_member,[nother 1]) - X_other;
        Xi_tenpermodel = repmat(X_single_member,[10*(nmodels-1) 1]) - X_tenpermodel_other;
        Xi_emean = squeeze(mean(reshape(Xi_tenpermodel,[nt 10*(nmodels-1) size(Xi,2)]),2));
        Xi_emean_f = lanczos_filter_datamatrix_trend_passthrough(Xi_emean,cutoff);
        
        Ct = cov(Xi_tenpermodel); 
        [tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs, pcvec, evl, rest] = forced_pattern_analysis(Xi, Xi_emean_f, truncation, scale,Ct);
        FPs       = insert_cols(FPs, icol_ret, icol_disc);
        
        tk_other(:,:,j) = tk(:,1:2);
        FPs_other(:,:,j) = FPs(1:2,:);
        s_other(:,j) = s;
        
        % construct data matrix for analysis of single model 
        npermodel = 30;
        Xi_single_model = [];
        members_other = [1:members(j)-1, members(j)+1:nei];
        members_other = members_other(randperm(nei-1));
        for k = 1:min([npermodel nei-1])
            ks = (members_other(k)-1)*nt+1:members_other(k)*nt;
            Xi_single_model = [Xi_single_model; X_single_member - X_single_model(ks,:)];
        end
        
        Ct = cov(Xi_single_model); 
        Xi_emean = squeeze(mean(reshape(Xi_single_model,[nt min([npermodel nei-1]) size(Xi,2)]),2));
        Xi_emean_f = lanczos_filter_datamatrix_trend_passthrough(Xi_emean,cutoff);
        [tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs, pcvec, evl, rest] = forced_pattern_analysis(Xi_single_model, Xi_emean_f, truncation, scale,Ct);
        FPs       = insert_cols(FPs, icol_ret, icol_disc);
        
        tk_single(:,:,j) = tk(:,1:2);
        FPs_single(:,:,j) = FPs(1:2,:);
        s_single(:,j) = s;
    end
    save([savedir,name,'_',num2str(cutoff),'yr_N',num2str(truncation),'_S2N_robustness_',models{i},'_other_models.mat'],'tk_other','FPs_other','s_other','members')
    save([savedir,name,'_',num2str(cutoff),'yr_N',num2str(truncation),'_S2N_robustness_',models{i},'_single_model.mat'],'tk_single','FPs_single','s_single','members')
    disp([models{i},' complete'])
end
