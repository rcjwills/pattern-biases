function [Xdisc_mean,Xdisc_anom,normalization] = plot_multivariate_model_bias_patterns(Xdisc,patternsf,ndisc,l,year1,lon,lat,scale1,scale2,ctrs1,ctrs2,title1,title2)

n = size(Xdisc,1);
ne = n/l;

dt = 1;

Xdisc_mean = zeros(l,1);

is_include = 1:ne;

t = (0:l-1)*dt+year1;

for i = is_include(end):-1:is_include(1)
    i1 = (i-1)*l+1;
    i2 = i*l;
    Xdisc_mean = Xdisc_mean + Xdisc(i1:i2,ndisc);
end
Xdisc_mean = Xdisc_mean./length(is_include);
normalization = std(Xdisc_mean);
Xdisc_mean = Xdisc_mean./normalization;
Xdisc = Xdisc./normalization;

if nargin > 5
    
    scale1 = scale1.*normalization;
    scale2 = scale2.*normalization;
    
    figure; subplot(8,1,[1 2 3]);
    try
        field = reshape(patternsf(ndisc,:),[2*length(lon) length(lat)]);
    catch
        lon = lon(1:length(lon)/2);
        field = reshape(patternsf(ndisc,:),[2*length(lon) length(lat)]);
    end
    field1 = scale1*field(1:end/2,:);
    field2 = scale2*field(end/2+1:end,:);
    plot_field_robinson_replace(lon,lat,field1,ctrs1);
    colorbar off;  hc1 = colorbar('eastoutside','position',[0.88    0.672    0.03    0.292]);
    %hc_position = get(hc1,'Position');
    %hc_position(4) = hc_position(4)*0.9;
    %set(hc1,'Position',hc_position)
    originalSize = get(gca, 'Position');
    set(gca,'fontsize',16)
    set(gcf,'position',[0 1000 500 700]);
    set(gca, 'Position', [originalSize(1)-0.11 originalSize(2) 1.25*originalSize(3) 1.25*originalSize(4)]);
    if nargin > 11
        title(title1,'fontsize',14)
    end
    subplot(8,1,[4 5 6]);
    plot_field_robinson_replace(lon,lat,field2,ctrs2);
    colorbar off;  hc2 = colorbar('eastoutside','position',[0.88    0.334    0.03    0.292]);
    %set(hc2,'Position',hc_position)
    originalSize = get(gca, 'Position');
    set(gca,'fontsize',16)
    set(gcf,'position',[0 1000 550 750]);
    set(gca, 'Position', [originalSize(1)-0.11 originalSize(2)-0.021 1.25*originalSize(3) 1.25*originalSize(4)]);
    if nargin > 11
        title(title2,'fontsize',14)
    end
    subplot(8,1,[7 8])
    
    for i = is_include(end):-1:is_include(1)
        i1 = (i-1)*l+1;
        i2 = i*l;
        hold on;
        plot(t,Xdisc(i1:i2,ndisc)-Xdisc_mean,'color',[0.7 0.7 0.7]);
        Xdisc_anom(:,i) = Xdisc(i1:i2,ndisc)-Xdisc_mean;
    end
    
    plot(t,Xdisc_mean,'k','linewidth',2);
    
    originalSize = get(gca, 'Position');
    originalSize(1) = originalSize(1)-0.04;
    set(gca, 'Position', originalSize);
    
    set(gca,'xlim',[floor(t(1)) ceil(t(end))]);
    set(gca,'ylim',[floor(min(Xdisc(:,ndisc))) ceil(max(Xdisc(:,ndisc)))])
    set(gca,'fontsize',14)
    xlabel('Year','fontsize',14);
    ylabel('Standard Deviations','fontsize',14);
    %xlabel('Year','fontsize',14);
else
    Xdisc_mean = Xdisc_mean./std(Xdisc_mean);
end