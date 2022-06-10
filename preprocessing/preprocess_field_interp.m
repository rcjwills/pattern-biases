function [data,clim,xtf,Ct,Mt,LON_AXIS,LAT_AXIS,AREA_WEIGHTS,icol_ret,icol_disc] = preprocess_field_interp(loni,lati,field,varnam,region,lon,lat)

s = size(field);

field(field>1e19) = nan;

%% Gridding and interpolation
sl = size(loni);
if sl(1)*sl(2)>(sl(1)+sl(2)) % ocean (2D) grid
    loni(loni<0) = loni(loni<0)+360;
    fieldi = field;
    field = zeros([length(lon) length(lat) s(3)]);
    [X,Y] = meshgrid(lon,lat);
    Xi = loni'; Yi = lati';
    %dispstat('','init')
    warning('off')
    if strcmp(varnam,'sic') || strcmp(varnam,'zos') || strcmp(varnam,'msftbarot') || strcmp(varnam,'ZOS_Monthly') || strcmp(varnam,'sos')
        disp('Interpolating . . .')
        for i = 1:s(3)
            %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
            tmp = griddata(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
            tmp2 = griddata(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y,'nearest')';
            tmp(isnan(tmp)) = tmp2(isnan(tmp)); % take nearest neighbor interpolation values when linear interpolation returns NaN's (near coast)
            field(:,:,i) = tmp;
        end
    else
        if sum(sum(sum(isnan(field)))) > 0
            disp('Warning: NaNs in data, interpolation will spread NaNs.')
        end
        disp('Interpolating . . .')
        for i = 1:s(3)
            %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
            field(:,:,i) = griddata(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
        end
    end
    warning('on')
    clear fieldi
elseif length(s) == 2 && length(loni) == length(lati) % unstructred grid
    fieldi = field;
    field = zeros([length(lon) length(lat) s(2)]);
    disp('Interpolating Unstructured Grid . . .')
    for i = 1:s(2)
        field(:,:,i) = griddata(loni,lati,fieldi(:,i),lon,lat')';
    end
    [X,Y] = meshgrid(lon,lat);
else % simple grid
    if length(loni)==length(lon) && length(lati)==length(lat)
        if sum(loni~=lon)==0 && sum(lati~=lat)==0
            disp('No interpolation needed')
            [X,Y] = meshgrid(lon,lat); % no interpolation needed
        else
            fieldi = field;
            field = zeros([length(lon) length(lat) s(3)]);
            [X,Y] = meshgrid(lon,lat);
            [Xi,Yi] = meshgrid(loni,lati);
            %dispstat('','init')
            warning('off')
            if strcmp(varnam,'sic') || strcmp(varnam,'zos')  || strcmp(varnam,'msftbarot') || strcmp(varnam,'ZOS_Monthly') || strcmp(varnam,'sos')
                disp('Interpolating . . .')
                for i = 1:s(3)
                    %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
                    tmp = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
                    tmp2 = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y,'nearest')';
                    tmp(isnan(tmp)) = tmp2(isnan(tmp)); % take nearest neighbor interpolation values when linear interpolation returns NaN's (near coast)
                    field(:,:,i) = tmp;
                end
            else
                if sum(sum(sum(isnan(field)))) > 0
                    disp('Warning: NaNs in data, interpolation will spread NaNs.')
                end
                disp('Interpolating . . .')
                for i = 1:s(3)
                    %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
                    field(:,:,i) = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
                end
            end
            warning('on')
            clear fieldi
        end
    else
        fieldi = field;
        field = zeros([length(lon) length(lat) s(3)]);
        [X,Y] = meshgrid(lon,lat);
        [Xi,Yi] = meshgrid(loni,lati);
        %dispstat('','init')
        warning('off')
        if strcmp(varnam,'sic') || strcmp(varnam,'zos')  || strcmp(varnam,'msftbarot') || strcmp(varnam,'ZOS_Monthly') || strcmp(varnam,'sos')
            disp('Interpolating . . .')
            for i = 1:s(3)
                %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
                tmp = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
                tmp2 = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y,'nearest')';
                tmp(isnan(tmp)) = tmp2(isnan(tmp)); % take nearest neighbor interpolation values when linear interpolation returns NaN's (near coast)
                field(:,:,i) = tmp;
            end
        else
            if sum(sum(sum(isnan(field)))) > 0
                disp('Warning: NaNs in data, interpolation will spread NaNs.')
            end
            disp('Interpolating . . .')
            for i = 1:s(3)
                %dispstat(sprintf(' . . . Interpolating . . . %d%%',round(i/s(3)*100)));
                field(:,:,i) = interp2(Xi,Yi,squeeze(fieldi(:,:,i))',X,Y)';
            end
        end
        warning('on')
        clear fieldi
    end
end
s = size(field);

tmp = mean(field,3);
%interpolate around lon = 360 or lon = 0 if there are NaNs
if isnan(sum(tmp(end,:)))
    for j = 1:s(2)
        for k = 1:s(3)
            field(end,j,k) = (field(end-1,j,k)+field(1,j,k))/2;
        end
    end
elseif isnan(sum(tmp(1,:)))
    for j = 1:s(2)
        for k = 1:s(3)
            field(1,j,k) = (field(end,j,k)+field(2,j,k))/2;
        end
    end
end
            
X = X'; Y = Y';
area = cos(Y*pi/180);

area(sum(isnan(field),3)./size(field,3)>0.4) = 0; % ignore grid point if more than 10% of data is missing

%% Remove ice area from surface temperatures

if strfind(varnam,'TS_sea')
    field(field<271.4) = 271.4;
end

%% compute anomaly from annual cycle
[anomalies,climatology] = monthly_anomalies(field);

%% domain masking

switch region
    case 'Pacific'
        domain = ones(size(area));
        domain(X<100) = 0;
        domain(X<103 & Y<5) = 0;
        domain(X<105 & Y<2) = 0;
        domain(X<111 & Y<-6) = 0;
        domain(X<114 & Y<-7) = 0;
        domain(X<127 & Y<-8) = 0;
        domain(X<147 & Y<-18) = 0;
        domain(Y>70) = 0;
        domain(Y>65 & (X<175 | X>200)) = 0;
        domain(Y<-45) = 0;
        domain(X>260 & Y>17) = 0;
        domain(X>270 & Y<=17 & Y>14) = 0;
        domain(X>276 & Y<=14 & Y>9) = 0;
        domain(X>290 & Y<=9) = 0;
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(icol_disc) = 0;
        
%     case 'Atlantic'
%         domain = ones(size(area));
%         domain(Y>70) = 0;
%         domain(Y<-45) = 0;
%         domain(X<260 & X>36 & Y>17) = 0;
%         domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
%         domain(X<276 & X>36 & Y<=14 & Y>9) = 0;
%         domain(X<290 & X>24 & Y<=9) = 0;
%         
%         domain = reshape(domain,s(1)*s(2),1)';
%         domain(icol_disc) = 0;
        
    case 'Atlantic_exp'
        domain = ones(size(area));
        domain(Y>85) = 0;
        domain(Y<-45) = 0;
        domain(X<260 & X>60 & Y>17) = 0;
        domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
        domain(X<276 & X>36 & Y<=14 & Y>9) = 0;
        domain(X<290 & X>24 & Y<=9) = 0;
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(icol_disc) = 0;
        
    case 'Atlantic'
        domain = ones(size(area));
        domain(Y>85) = 0;
        domain(Y<0) = 0;
        domain(X<260 & X>60 & Y>17) = 0;
        domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
        domain(X<276 & X>36 & Y<=14 & Y>9) = 0;
        domain(X<290 & X>24 & Y<=9) = 0;
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(icol_disc) = 0;
        
%     case 'Atlantic' 
%         domain = ones(size(area));
%         domain(Y>85) = 0;
%         domain(Y<-55) = 0;
%         domain(X<260 & X>100 & Y>17) = 0;
%         domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
%         domain(X<276 & X>36 & Y<=14 & Y>9) = 0;
%         domain(X<290 & X>24 & Y<=9) = 0;
%         domain(X<180 & Y>10 & Y<50) = 0; % exclude Med
%         domain(X>355 & Y>30 & Y<40) = 0; % exclude Med
%         domain(X<290 & X>100 & Y<70 & Y>50) = 0; % exclude Hudson 
%         domain(X<279 & X>100 & Y>70) = 0;
%         domain(X<330 & X>100 & Y>80) = 0;
%         domain(X<100 & X>9 & Y>50 & Y<=56) = 0; % exclude Baltic
%         domain(X<100 & X>15 & Y>56 & Y<=66) = 0; % exclude Baltic
%         
%         domain = reshape(domain,s(1)*s(2),1)';
%         domain(icol_disc) = 0;
        
    case 'Pacific_20-70N'
        domain = ones(size(area));
        domain(X<105) = 0;
        domain(Y>70) = 0;
        domain(Y<20) = 0;
        domain(X>260) = 0;
        domain(Y>65 & (X<175 | X>200)) = 0;
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(icol_disc) = 0;
        
    case 'Global'
        domain = ones(size(area));
        domain = reshape(domain,s(1)*s(2),1)';
        
    case 'Global_Ocean'
        domain = ones(size(area));
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(icol_disc) = 0;
        
    case 'Southern_Ocean'
        domain = ones(size(area));
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(Y>-35) = 0;
        domain(icol_disc) = 0;
        
    case 'AMOC'
        domain = ones(size(area));
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(X<-35) = 0;
        domain(icol_disc) = 0;
    case 'PMOC'
        domain = ones(size(area));
        
        domain = reshape(domain,s(1)*s(2),1)';
        domain(X<-35) = 0;
        domain(icol_disc) = 0;
end
    
data = anomalies;

s = size(data);
x = reshape(data,s(1)*s(2),s(3))';
AREA_WEIGHTS = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';
Mt = reshape(climatology,s(1)*s(2),12)';
LON_AXIS = lon;
LAT_AXIS = lat;

icol_ret = find(AREA_WEIGHTS~=0 & domain);
icol_disc = find(AREA_WEIGHTS==0 | ~domain);
xtf = x(:,icol_ret);
AREA_WEIGHTS = AREA_WEIGHTS(icol_ret);
Mt = Mt(:,icol_ret);

data = reshape(insert_cols(xtf,icol_ret,icol_disc)',[s(1) s(2) s(3)]);
clim = reshape(insert_cols(Mt,icol_ret,icol_disc)',[s(1) s(2) 12]);

%disp('Computing covariance matrix...')
%Ct               = cov(xtf);
Ct = nan;
