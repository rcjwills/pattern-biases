function [] = preprocess_cmip5_LE(model_name,year1,year2,averaging,realm1,var1,realm2,var2)

% load two LE variables, compute local trends and leading EOFs

%% Input Ensemble and Variable Information

% EXAMPLE INPUT
% model_name = 'cesm1';
% year1 = 1979;
% year2 = 2020;
% averaging = 'Amean';
% 
% 
% realm1 = 'Omon';
% var1 = 'tos';
% realm2 = 'Amon';
% var2 = 'psl';

% flag for turning ts into sst
if strcmp(var1,'ts-sst')
    if strcmp(var2,'ts')
        disp('Error: Cannot have ts-sst for both variables')
    end
    ts_to_sst = 1;
    var1 = 'ts';
elseif strcmp(var2,'ts-sst')
    ts_to_sst = 1;
    var2 = 'ts';
else
    ts_to_sst = 0;
end

LON_AXIS = (1.25:2.5:358.75)';
LAT_AXIS = (-88.75:2.5:88.75)';

dir = ['/eos3/rcwills/cmip_ensembles/',model_name,'_lens/'];

save_dir = '/eos15/rcwills/projects_output/MMLEA/';
save_name_trends = [model_name,'_',var1,'_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
save_name_datamatrix = [model_name,'_',var1,'_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
if ts_to_sst
    if strcmp(var1,'ts')
        save_name_trends = [model_name,'_ts-sst_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
        save_name_datamatrix = [model_name,'_ts-sst_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
    elseif strcmp(var2,'ts')
        save_name_trends = [model_name,'_',var1,'_ts-sst_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
        save_name_datamatrix = [model_name,'_',var1,'_ts-sst_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
    end
end    

%% Determine Ensemble Size and File Structure

files1 = get_files([dir,realm1,'/',var1,'/']);
files1 = select_files(files1,'historical_rcp85');
files2 = get_files([dir,realm2,'/',var2,'/']);
files2 = select_files(files2,'historical_rcp85');

% modify based on lon/lat variable names
try
    [lon1,lat1,field1] = get_avg_field_3d(files1{1},{'lon','lat',var1});
catch
    try
        [lon1,lat1,field1] = get_avg_field_3d(files1{1},{'longitude','latitude',var1});
    catch
        [lon1,lat1,field1] = get_avg_field_3d(files1{1},{'xt_ocean','yt_ocean',var1});
    end
end
[lon2,lat2,field2] = get_avg_field_3d(files2{1},{'lon','lat',var2});

if length(files2)~=length(files1)
    disp('Error: different number of files for each variable')
end

member = zeros(length(files1),1);

for i = 1:length(files1)
    k1 = strfind(files1{i},'_r'); k1 = k1(end);
    k2 = strfind(files1{i},'i1');
    member(i) = str2double(files1{i}(k1+2:k2-1));
end

nfiles = length(find(member==1));
ne = length(files1)./nfiles;

y1(nfiles) = 0;
y2(nfiles) = 0;

for j = 1:nfiles
    file = files1{j+nfiles}; % using second ensemble member because some models have differences in first ensemble member
    k3 = strfind(file,'p1_');
    y1(j) = str2double(file(k3+3:k3+6));
    y2(j) = str2double(file(k3+10:k3+13));
end

% find which files are needed for time period of interest
js = find(y1<=year1,1,'last'):find(y2>=year2,1,'first');
y0 = y1(js(1));
yend = y2(js(end));

nt = 12*(year2-year1+1);
ncut0 = 12*(year1-y0);
ncutend = 12*(yend-year2);

if strcmp(model_name,'cesm1') % models where first ensemble member is a different length
    for j = 1:nfiles
        file = files1{j}; % first ensemble member
        k3 = strfind(file,'p1_');
        y1(j) = str2double(file(k3+3:k3+6));
        y2(j) = str2double(file(k3+10:k3+13));
    end
    
    % find which files are needed for time period of interest
    js = find(y1<=year1,1,'last'):find(y2>=year2,1,'first');
    y0 = y1(js(1));
    yend = y2(js(end));
    
    ncut0_1 = 12*(year1-y0);
    ncutend_1 = 12*(yend-year2);
else
    ncut0_1 = ncut0;
    ncutend_1 = ncutend;
end

%% 

xtf1_all = [];
xtf2_all = [];
Mt1_emean = 0;
Mt2_emean = 0;

field1_trends = zeros([length(LON_AXIS),length(LAT_AXIS),ne]);
field2_trends = zeros([length(LON_AXIS),length(LAT_AXIS),ne]);

for i = 1:ne
    disp(['Processing Ensemble Member ',num2str(i)])
    is = (i-1)*nfiles+1:i*nfiles;
    is = is(js);
    
    if i == 1
        field1 = concatenate_fields(files1(is),var1,0);
        field1 = field1(:,:,ncut0_1+1:end-ncutend_1);
        field2 = concatenate_fields(files2(is),var2,0);
        field2 = field2(:,:,ncut0_1+1:end-ncutend_1);
    else
        field1 = concatenate_fields(files1(is),var1,0);
        field1 = field1(:,:,ncut0+1:end-ncutend);
        field2 = concatenate_fields(files2(is),var2,0);
        field2 = field2(:,:,ncut0+1:end-ncutend);
    end
    
    if strcmp(var1,'tos')
        field1(abs(field1)>1e10) = nan;
    end
    
    if ts_to_sst
        if strcmp(var1,'ts')
            field1(field1<271.4) = 271.4;
            if i == 1
                [Y,X] = meshgrid(lat1,lon1);
                X(X>180) = X(X>180)-360;
                land = landmask(Y,X);
                land = repmat(land,[1 1 size(field1,3)]);
            end
            field1(land) = nan;
        elseif strcmp(var2,'ts')
            field2(field2<271.4) = 271.4;
            if i == 1
                [Y,X] = meshgrid(lat2,lon2);
                X(X>180) = X(X>180)-360;
                land = landmask(Y,X);
                land = repmat(land,[1 1 size(field2,3)]);
            end
            field2(land) = nan;
        end
    end
            
    months = repmat(1:12,[1 nt/12]);
    years = reshape(repmat(year1:year2,[12 1]),[nt 1])';
    t = year1:year2; t= t';
    switch averaging
        case 'Amean'
            field1 = monthly_to_seasonal(months,years,field1,1:12,3);
            field2 = monthly_to_seasonal(months,years,field2,1:12,3);
        case 'DJF'
            field1 = monthly_to_seasonal(months,years,field1,[12 1 2],3);
            field2 = monthly_to_seasonal(months,years,field2,[12 1 2],3);
            t = t(2:end);
        case 'JJA'
            field1 = monthly_to_seasonal(months,years,field1,6:8,3);
            field2 = monthly_to_seasonal(months,years,field2,6:8,3);
        otherwise
            t = year1:1/12:year2+11/12; t = t'; % monthly
    end
    
    %% interpolate data and compute covariance matrix, field 1
    % note, the local version of preprocess_field_interp does not calculate
    % the covariance matrix
    [data,~,xtf1,~,Mt1,LON_AXIS,LAT_AXIS,AREA_WEIGHTS1,icol_ret1,icol_disc1] = preprocess_field_interp(lon1,lat1,field1,var1,'Global',LON_AXIS,LAT_AXIS);
    
    [~,field1_trends(:,:,i)] = detrend(data,3,1,t);
    
    %% interpolate data and compute covariance matrix, field 2
    [data,~,xtf2,~,Mt2,LON_AXIS,LAT_AXIS,AREA_WEIGHTS2,icol_ret2,icol_disc2] = preprocess_field_interp(lon2,lat2,field2,var2,'Global',LON_AXIS,LAT_AXIS);
    
    [~,field2_trends(:,:,i)] = detrend(data,3,1,t);
    
    Mt1_emean = Mt1_emean+Mt1./ne;
    Mt2_emean = Mt2_emean+Mt2./ne;
    xtf1_all = [xtf1_all; xtf1];
    xtf2_all = [xtf2_all; xtf2];
end

%%

save([save_dir,save_name_trends],'field1_trends','field2_trends','lon1','lat1','lon2','lat2')
save([save_dir,save_name_datamatrix],'xtf1_all','xtf2_all','Mt1_emean','Mt2_emean','LON_AXIS','LAT_AXIS','icol_ret1','icol_disc1','AREA_WEIGHTS1','icol_ret2','icol_disc2','AREA_WEIGHTS2','-v7.3')