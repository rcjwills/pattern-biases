function [h,hc,landareas] = plot_field_robinson(lon,lat,field,ctrs,origin,latlimit,lonlimit,wrap,filter,cmap)

field = squeeze(field);

%colors = {[0 0 1] [0 0.5 1] [1 0.5 0] [0.75 0 0]}; %blue-red
%colors = {[0.6 0 0] [1 0.5 0] [0 0.6 1] [0 0 0.6]}; %red-blue
%colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [49 123 184]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [189 34 39]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[183 106 41]/256 [225 165 100]/256 [245 224 158]/256 [164 213 169]/256 [110 170 200]/256 [7 57 87]/256}; %precip

if nargin > 8
    try
        field = conv2(field,filter,'same');
    catch
        disp('No filter applied')
    end
end

if nargin > 7
    if strcmp(wrap,'wrap')
        field = [field; field(1,:)];
        try
            lon = [lon; lon(1)+360];
        catch
            lon = [lon, lon(1)+360];
        end
    end
end

if nargin < 4
    ctrs = linspace(min(min(field)),max(max(field)),20);
end

if nargin < 7
    origin = [0 180 0];
    latlimit = [-90 90];
    lonlimit = [0 360];
end

field(field<min(ctrs)) = min(ctrs);

lcmap = length(ctrs)-1;
ch = [min(ctrs)-1e-8 max(ctrs)];

figure;
axesm ('robinson', 'Frame', 'on', 'Grid', 'off','origin',origin,'maplatlimit',latlimit,'maplonlimit',lonlimit);
[~,h] = contourfm(lat,lon,field',ctrs,'linestyle','none');
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'EdgeColor',[0.15 0.15 0.15],'facecolor',[1 1 1],'facealpha',0.15);
caxis(ch); hc = colorbar; 
if nargin > 9
    colormap(cmap);
else
    cmap = colormap(diverging0(colors,[1 1 1],lcmap));
end
pretty_figure(700,300,'none','none','none','none',16);  axis off
set(gcf,'renderer','zbuffer')

% edge_darkening = 0.7;
% hold on;
% if min(ctrs) == -max(ctrs)
%     for i = 1:length(cmap)/2
%         hold on; contour(lon,lat',field',[ctrs(i) ctrs(i)],'LineColor',cmap(i,:)*edge_darkening);
%     end
%     for i = length(cmap)/2+1:length(cmap)
%         if sum(cmap(i,:)==1)~=3
%             hold on; contour(lon,lat',field',[ctrs(i) ctrs(i)],'LineColor',cmap(i,:)*edge_darkening);
%         end
%     end
% else
%     for i = 1:length(cmap)
%         if sum(cmap(i,:)==1)~=3
%             hold on; contour(lon,lat',field',[ctrs(i) ctrs(i)],'LineColor',cmap(i,:)*edge_darkening);
%         end
%     end
% end