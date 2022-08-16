function [rescaling] = trend_comparison_figures(lon,lat,trend_lens,trend_obs,ctrs,ctrs_std,rescale,model_name,obs_name,label,trend_std_input)

save_flag = 0;

if max(ctrs)>20
    save_dir = '/Specify Save Directory for SLP Figures/';
else
    save_dir = '/Specify Save Directory for SST Figures/';
end

if rescale==2 % SST, rescale global mean of models to match global mean of observations
    trend_std = std(trend_lens,0,3);
    trend_lens_rescaled = mean(trend_lens,3)./global_mean(lon,lat,mean(trend_lens,3)).*global_mean(lon,lat,trend_obs);
    rescaling = global_mean(lon,lat,trend_obs)./global_mean(lon,lat,mean(trend_lens,3));
elseif rescale==1 % SLP, remove global mean
    for i = 1:size(trend_lens,3)
        trend_lens(:,:,i) = trend_lens(:,:,i) - global_mean(lon,lat,trend_lens(:,:,i));
    end
    trend_lens_rescaled = mean(trend_lens,3);
    trend_std = std(trend_lens,0,3);
    trend_obs = trend_obs - global_mean(lon,lat,trend_obs);
else
    trend_std = std(trend_lens,0,3);
    trend_lens_rescaled = mean(trend_lens,3);
end

if nargin > 10
    trend_std = trend_std_input;
end

origin = [0 180 0];
latlim = [-80 80];
lonlim = [0 360];

% origin = [0 -90 0];
% latlim = [-80 80];
% lonlim = [-270 90];

if(1)
    plot_field_robinson(lon,lat,trend_obs,ctrs,origin,latlim,lonlim,'wrap');
    switch max(ctrs)
        case 1.6
            hc = colorbar; set(hc,'xtick',-1.6:0.8:1.6); set(hc,'xticklabel',{'-1.6','-0.8','0','0.8','1.6°C'})
        case 240
            hc = colorbar; set(hc,'xtick',-200:100:200); set(hc,'xticklabel',{'-200','-100','0','100','200 Pa'})
        case 400
            hc = colorbar; set(hc,'xtick',-400:200:400); set(hc,'xticklabel',{'-400','-200','0','200','400 Pa'})
        case 600
            hc = colorbar; set(hc,'xtick',-600:200:600); set(hc,'xticklabel',{'-600','-400','-200','0','200','400','600 Pa'})
    end
    if save_flag==1
        saveas(gcf,[save_dir,obs_name,'_',label,'.eps'],'epsc')
    end
end

if nargin > 10
    if(0)
        plot_field_robinson(lon,lat,trend_std_input,ctrs,origin,latlim,lonlim,'wrap');
        switch max(ctrs)
            case 1.6
                hc = colorbar; set(hc,'xtick',-1.6:0.8:1.6); set(hc,'xticklabel',{'-1.6','-0.8','0','0.8','1.6°C'})
            case 240
                hc = colorbar; set(hc,'xtick',-200:100:200); set(hc,'xticklabel',{'-200','-100','0','100','200 Pa'})
            case 400
                hc = colorbar; set(hc,'xtick',-400:200:400); set(hc,'xticklabel',{'-400','-200','0','200','400 Pa'})
            case 600
                hc = colorbar; set(hc,'xtick',-600:200:600); set(hc,'xticklabel',{'-600','-400','-200','0','200','400','600 Pa'})
        end
    end
end

plot_field_robinson(lon,lat,trend_lens_rescaled,ctrs,origin,latlim,lonlim,'wrap');
switch max(ctrs)
    case 1.6
        hc = colorbar; set(hc,'xtick',-1.6:0.8:1.6); set(hc,'xticklabel',{'-1.6','-0.8','0','0.8','1.6°C'})
    case 240
        hc = colorbar; set(hc,'xtick',-200:100:200); set(hc,'xticklabel',{'-200','-100','0','100','200 Pa'})
    case 400
        hc = colorbar; set(hc,'xtick',-400:200:400); set(hc,'xticklabel',{'-400','-200','0','200','400 Pa'})
    case 600
        hc = colorbar; set(hc,'xtick',-600:200:600); set(hc,'xticklabel',{'-600','-400','-200','0','200','400','600 Pa'})
end
title(model_name)
if save_flag==1
    saveas(gcf,[save_dir,model_name,'_',label,'.eps'],'epsc')
end

field = (trend_obs-trend_lens_rescaled)./trend_std;
field(isinf(field)) = nan;
[Y,X] = meshgrid(lat,lon);

plot_field_robinson(lon,lat,field,ctrs_std,origin,latlim,lonlim,'wrap');
hc = colorbar; set(hc,'xtick',-6:2:6);
hold on; contourm(lat,lon,field',[2 2],'k')
hold on; contourm(lat,lon,field',[-2 -2],'k')
title(model_name)
annotation('textbox',[0.27 0.83 .1 .1],'String',['N = ',num2str(size(trend_lens,3))],'fontsize',14,'EdgeColor','none')
field(abs(Y)>=60) = nan; % only compute RMSE over [60S 60N]
annotation('textbox',[0.61 0.83 .25 .1],'String',['RMSE = ',num2str(round(sqrt(global_mean(lon,lat,field.^2)),2)),' \sigma'],'fontsize',14,'EdgeColor','none')
if save_flag==1
    saveas(gcf,[save_dir,'stdev_',model_name,'_',obs_name,'_',label,'.eps'],'epsc')
end

