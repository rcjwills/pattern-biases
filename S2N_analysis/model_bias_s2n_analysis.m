function [tk, FPs, LON_AXIS, LAT_AXIS, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs] = model_bias_s2n_analysis(LON_AXIS,LAT_AXIS,field_lens,field_obs,truncation,cutoff,time_option)

n = size(field_lens,3); 
nyr = n/12;
months = repmat(1:12,[1 nyr]);
years = ceil(1/12:1/12:nyr);
ne = size(field_lens,4);

[Y,X] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
switch time_option
    case 'monthly'
        [X2,Y2,~,field_obs] = coarsen2(X,Y,area,field_obs);
        for i = 1:ne
            tmp = field_lens(:,:,:,i);
            [~,~,~,field_lens_new(:,:,:,i)] = coarsen2(X,Y,area,tmp);
        end
    case 'Amean'
        tmp = monthly_to_seasonal(months,years,field_obs,1:12,3);
        [X2,Y2,~,field_obs] = coarsen2(X,Y,area,tmp);
        for i = 1:ne
            tmp = monthly_to_seasonal(months,years,field_lens(:,:,:,i),1:12,3);
            [~,~,~,field_lens_new(:,:,:,i)] = coarsen2(X,Y,area,tmp);
        end
        n = nyr;
    case 'DJF'
        tmp = monthly_to_seasonal(months,years,field_obs,[12 1 2],3);
        [X2,Y2,~,field_obs] = coarsen2(X,Y,area,tmp);
        for i = 1:ne
            tmp = monthly_to_seasonal(months,years,field_lens(:,:,:,i),[12 1 2],3);
            [~,~,~,field_lens_new(:,:,:,i)] = coarsen2(X,Y,area,tmp);
        end
        n = nyr-1;
    case 'JJA'
        tmp = monthly_to_seasonal(months,years,field_obs,6:8,3);
        [X2,Y2,~,field_obs] = coarsen2(X,Y,area,tmp);
        for i = 1:ne
            tmp = monthly_to_seasonal(months,years,field_lens(:,:,:,i),6:8,3);
            [~,~,~,field_lens_new(:,:,:,i)] = coarsen2(X,Y,area,tmp);
        end
        n = nyr;
end      
field_lens = field_lens_new; clear tmp field_lens_new
LON_AXIS = X2(:,1); clear X2
LAT_AXIS = Y2(1,:); clear Y2

nlon = size(field_lens,1);
nlat = size(field_lens,2);
t = 1:n;

%% Parameters
% truncation = number of EOFs retained (e.g., 10-50)
% cutoff - focus on variability at timescales longer than this cutoff (in
% timesteps, e.g., 60 months)

%% Preprocessing

% remove time mean
field_lens = field_lens - repmat(mean(mean(field_lens,4),3),[1 1 n ne]);
field_obs = field_obs - repmat(mean(field_obs,3),[1 1 n]);

% remove ensemble mean and reshape data
lens_emean = mean(field_lens,4);
lens_anomalies = field_lens - repmat(lens_emean,[1 1 1 ne]);
obs_anomalies = field_obs - lens_emean;
anomalies_all = lens_anomalies;
anomalies_all(:,:,:,ne+1) = obs_anomalies;
anomalies_all = reshape(anomalies_all,[nlon nlat n*(ne+1)]); % concatenate ensemble members in time

%s = size(sst_anomalies_all);
% for i = 1:s(1)
%     for j = 1:s(2)
%         for k = 1:s(4)
%             p = polyfit(time,squeeze(sst_anomalies_all(i,j,:,k))',1);
%             sst_trends_all(i,j,k) = p(1);
%         end
%         p = polyfit(time,squeeze(ERSST_anomalies(i,j,:))',1);
%         ERSST_trend(i,j) = p(1);
%     end
% end

[Y,~] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
area(isnan(mean(anomalies_all,3))) = 0; 

domain = ones(size(area));

% Can specify sub-domains to analyze using the following form:

%domain(Y>60) = 0;
%domain(Y<-60) = 0;

s = size(anomalies_all);
X = reshape(anomalies_all,s(1)*s(2),s(3))'; 
Xe = reshape(obs_anomalies,s(1)*s(2),n)';
%s = size(sst_trends_all);
%X = reshape(sst_trends_all,s(1)*s(2),s(3))';
%Xe = reshape(ERSST_trend,s(1)*s(2),1)';
AREA_WEIGHTS = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret = find(AREA_WEIGHTS~=0 & domain);
icol_disc = find(AREA_WEIGHTS==0 | ~domain);
X = X(:,icol_ret);
Xe = Xe(:,icol_ret);
AREA_WEIGHTS = AREA_WEIGHTS(icol_ret);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS' ./ sum(AREA_WEIGHTS);
scale    = sqrt(normvec);

Xe_f = Xe; t = t';
for i = 1:size(Xe,2)
    p = polyfit(t,Xe(:,i),1);
    tmp = Xe(:,i)-p(1)*t-p(2);
    tmp = lanczos([flipud(tmp); tmp; flipud(tmp)],1,cutoff);
    Xe_f(:,i) = tmp((end/3+1):2*end/3)+p(1)*t+p(2);
end

%% Calculate forced patterns / Signal-to-noise maximizing EOF analysis
[tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs] = forced_pattern_analysis(X, Xe_f, truncation, scale);
FPs       = insert_cols(FPs, icol_ret, icol_disc);
EOFs       = insert_cols(EOFs, icol_ret, icol_disc);
fingerprintsf        = insert_rows(fingerprints, icol_ret, icol_disc);

