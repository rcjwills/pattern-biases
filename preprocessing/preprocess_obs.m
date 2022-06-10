function [trend,LAT_AXIS,LON_AXIS,xtf,icol_ret,icol_disc] = preprocess_obs(lon,lat,field,year0,year1,year2,averaging)

% year0 is the first year of the dataset, year1 and year2 specify the range
% to analyze

LON_AXIS = (1.25:2.5:358.75)';
LAT_AXIS = (-88.75:2.5:88.75)';

yend = year0+size(field,3)./12-1;

nt = 12*(year2-year1+1);
ncut0 = 12*(year1-year0);
ncutend = 12*(yend-year2);

field = field(:,:,ncut0+1:end-ncutend);

% [Y,X] = meshgrid(lat,lon);
% X(X>180) = X(X>180)-360;
% land = landmask(Y,X);
% land = repmat(land,[1 1 size(field,3)]);
% field(land) = nan;

months = repmat(1:12,[1 nt/12]);
years = reshape(repmat(year1:year2,[12 1]),[nt 1])';
t = year1:year2; t= t';
switch averaging
    case 'Amean'
        field = monthly_to_seasonal(months,years,field,1:12,3);
    case 'DJF'
        field = monthly_to_seasonal(months,years,field,[12 1 2],3);
        t = t(2:end);
    case 'JJA'
        field = monthly_to_seasonal(months,years,field,6:8,3);
    otherwise
        t = year1:1/12:year2+11/12; t = t'; % monthly
end

[data,~,xtf,~,Mt,LON_AXIS,LAT_AXIS,AREA_WEIGHTS,icol_ret,icol_disc] = preprocess_field_interp(lon,lat,field,'varnam','Global',LON_AXIS,LAT_AXIS);
    
[~,trend] = detrend(data,3,1,t);