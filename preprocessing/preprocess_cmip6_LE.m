function [] = preprocess_cmip6_LE(model_name,ssp_name,members,year1,year2,averaging,realm1,var1,realm2,var2)

% load two LE variables, interpolate if needed, compute local trends and
% save datamatrix

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
ts_to_sst = 0;

LON_AXIS = (1.25:2.5:358.75)';
LAT_AXIS = (-88.75:2.5:88.75)';

dir = ['/eos3/rcwills/cmip_ensembles/',model_name,'_lens/'];

save_dir = '/eos15/rcwills/projects_output/MMLEA/';
save_name_trends = [model_name,'_',var1,'_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
save_name_datamatrix = [model_name,'_',var1,'_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
if ts_to_sst
    if strcmp(var1,'ts')
        if strcmp(var2,'ts')
            disp('Error: Cannot have ts for both variables')
        end
        save_name_trends = [model_name,'_ts-sst_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
        save_name_datamatrix = [model_name,'_ts-sst_',var2,'_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
    elseif strcmp(var2,'ts')
        save_name_trends = [model_name,'_',var1,'_ts-sst_',num2str(year1),'-',num2str(year2),'_',averaging,'_trends.mat'];
        save_name_datamatrix = [model_name,'_',var1,'_ts-sst_',num2str(year1),'-',num2str(year2),'_',averaging,'_datamatrix.mat'];
    end
end    

if strcmp(model_name,'cesm2')
    if strcmp(var1,'tos')
        var1 = 'SST';
    end
    if strcmp(var2,'psl')
        var2 = 'PSL';
    end
    hist_name = 'HIST';
    ssp_name = 'SSP370';
else
    hist_name = 'historical_';
end

%% Determine Ensemble Size and File Structure

files1 = get_files([dir,realm1,'/',var1,'/']);
files1_hist = select_files(files1,hist_name);
files1_ssp = select_files(files1,ssp_name);
files2 = get_files([dir,realm2,'/',var2,'/']);
files2_hist = select_files(files2,hist_name);
files2_ssp = select_files(files2,ssp_name);

if strcmp(model_name,'canesm5')
    files1_hist = select_files(files1_hist,'p2f1');
    files1_ssp = select_files(files1_ssp,'p2f1');
    files2_hist = select_files(files2_hist,'p2f1');
    files2_ssp = select_files(files2_ssp,'p2f1');
end

% modify based on lon/lat variable names
try
    [lon1,lat1,field1] = get_avg_field_3d(files1_hist{1},{'lon','lat',var1});
catch
    try
        [lon1,lat1,field1] = get_avg_field_3d(files1_hist{1},{'longitude','latitude',var1});
    catch
        [lon1,lat1,field1] = get_avg_field_3d(files1_hist{1},{'xt_ocean','yt_ocean',var1});
    end
end
[lon2,lat2,field2] = get_avg_field_3d(files2_hist{1},{'lon','lat',var2});

% if length(files2_hist)~=length(files1_hist) || length(files2_ssp)~=length(files1_ssp)
%     disp('Warning: different number of files for each variable')
% end

member_hist_1 = zeros(length(files1_hist),1);
member_ssp_1 = zeros(length(files1_ssp),1);
member_hist_2 = zeros(length(files1_hist),1);
member_ssp_2 = zeros(length(files1_ssp),1);

if strcmp(model_name,'giss_e21g')
    for i = 1:length(files1_hist)
        k1 = strfind(files1_hist{i},'_r'); k1 = k1(end);
        k2 = strfind(files1_hist{i},'i1');
        member_hist_1(i) = str2double([files1_hist{i}(k1+2:k2-1),files1_hist{i}(k2+3),files1_hist{i}(k2+5)]);
    end
    
    for i = 1:length(files1_ssp)
        k1 = strfind(files1_ssp{i},'_r'); k1 = k1(end);
        k2 = strfind(files1_ssp{i},'i1');
        member_ssp_1(i) = str2double([files1_ssp{i}(k1+2:k2-1),files1_ssp{i}(k2+3),files1_ssp{i}(k2+5)]);
    end
    
    for i = 1:length(files2_hist)
        k1 = strfind(files2_hist{i},'_r'); k1 = k1(end);
        k2 = strfind(files2_hist{i},'i1');
        member_hist_2(i) = str2double([files2_hist{i}(k1+2:k2-1),files2_hist{i}(k2+3),files2_hist{i}(k2+5)]);
    end
    
    for i = 1:length(files2_ssp)
        k1 = strfind(files2_ssp{i},'_r'); k1 = k1(end);
        k2 = strfind(files2_ssp{i},'i1');
        member_ssp_2(i) = str2double([files2_ssp{i}(k1+2:k2-1),files2_ssp{i}(k2+3),files2_ssp{i}(k2+5)]);
    end
elseif strcmp(model_name,'cesm2')
    for i = 1:length(files1_hist)
        k1 = strfind(files1_hist{i},'LE2');
        member_hist_1(i) = str2double([files1_hist{i}(k1+4:k1+7),files1_hist{i}(k1+10:k1+11)]);
    end
    
    for i = 1:length(files1_ssp)
        k1 = strfind(files1_ssp{i},'LE2');
        member_ssp_1(i) = str2double([files1_ssp{i}(k1+4:k1+7),files1_ssp{i}(k1+10:k1+11)]);
    end
    
    for i = 1:length(files2_hist)
        k1 = strfind(files2_hist{i},'LE2');
        member_hist_2(i) = str2double([files2_hist{i}(k1+4:k1+7),files2_hist{i}(k1+10:k1+11)]);
    end
    
    for i = 1:length(files2_ssp)
        k1 = strfind(files2_ssp{i},'LE2');
        member_ssp_2(i) = str2double([files2_ssp{i}(k1+4:k1+7),files2_ssp{i}(k1+10:k1+11)]);
    end
else
    for i = 1:length(files1_hist)
        k1 = strfind(files1_hist{i},'_r'); k1 = k1(end);
        k2 = strfind(files1_hist{i},'i1');
        member_hist_1(i) = str2double(files1_hist{i}(k1+2:k2-1));
    end
    
    for i = 1:length(files1_ssp)
        k1 = strfind(files1_ssp{i},'_r'); k1 = k1(end);
        k2 = strfind(files1_ssp{i},'i1');
        member_ssp_1(i) = str2double(files1_ssp{i}(k1+2:k2-1));
    end
    
    for i = 1:length(files2_hist)
        k1 = strfind(files2_hist{i},'_r'); k1 = k1(end);
        k2 = strfind(files2_hist{i},'i1');
        member_hist_2(i) = str2double(files2_hist{i}(k1+2:k2-1));
    end
    
    for i = 1:length(files2_ssp)
        k1 = strfind(files2_ssp{i},'_r'); k1 = k1(end);
        k2 = strfind(files2_ssp{i},'i1');
        member_ssp_2(i) = str2double(files2_ssp{i}(k1+2:k2-1));
    end
end


nfiles_hist = length(find(member_hist_2==member_hist_2(1)));
nfiles_ssp = length(find(member_ssp_2==member_ssp_2(1)));
ne = length(members);

nfiles = nfiles_hist+nfiles_ssp;

clear files1 files2

for i = 1:ne
    is_hist = (i-1)*nfiles+1:(i-1)*nfiles+nfiles_hist;
    is_ssp = (i-1)*nfiles+nfiles_hist+1:i*nfiles;
    files1{is_hist} = files1_hist{member_hist_1==members(i)};
    files1{is_ssp} = files1_ssp{member_ssp_1==members(i)};
    files2{is_hist} = files2_hist{member_hist_2==members(i)};
    files2{is_ssp} = files2_ssp{member_ssp_2==members(i)};
end

y1(nfiles) = 0;
y2(nfiles) = 0;

% Note, assuming no gap between historical and ssp

js = 1:nfiles; % including all files until there is a need to do otherwise
% needs to be modified if time variable is not days since 01.01.1850
time = get_avg_field_3d(files1_hist{nfiles_hist+1},'time');
if strcmp(model_name,'norcpm1')
    y0 = floor(time(1)./365)+1800; % days since 1800 for some reason
else
    y0 = floor(time(1)./365)+1850;
end
time = get_avg_field_3d(files1_ssp{end},'time');
if strcmp(model_name,'cesm2')
    yend = floor(time(end)./365)+2015-1;
elseif strcmp(model_name,'norcpm1')
    yend = 2029;
else
    yend = floor(time(end)./365)+1850-1;
end

nt = 12*(year2-year1+1);
ncut0 = 12*(year1-y0);
ncutend = 12*(yend-year2);

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

    if strcmp(model_name,'cnrm_cm6') % has different lengths of simulations
        time = get_avg_field_3d(files1{is(1)},'time');
        y0 = floor(time(1)./365)+1850;
        time = get_avg_field_3d(files1{is(end)},'time');
        yend = floor(time(end)./365)+1850-1;
        ncut0 = 12*(year1-y0);
        ncutend = 12*(yend-year2);
    end
 
    field1 = concatenate_fields(files1(is),var1,0);
    field1 = field1(:,:,ncut0+1:end-ncutend);
    field2 = concatenate_fields(files2(is),var2,0);
    field2 = field2(:,:,ncut0+1:end-ncutend);
    
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

save([save_dir,save_name_trends],'field1_trends','field2_trends','lon1','lat1','lon2','lat2','ssp_name','members')
save([save_dir,save_name_datamatrix],'xtf1_all','xtf2_all','Mt1_emean','Mt2_emean','LON_AXIS','LAT_AXIS','icol_ret1','icol_disc1','AREA_WEIGHTS1','icol_ret2','icol_disc2','AREA_WEIGHTS2','ssp_name','members','-v7.3')
