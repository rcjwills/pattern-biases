function [Xdisc_mean] = plot_ensemble_and_obs_patterns(Xdisc,patternsf,ndisc,l,year1,lon,lat,ctrs,title_text)

if nargin > 5
    figure; subplot(3,1,[1 2]);
    field = reshape(patternsf(ndisc,:),[length(lon) length(lat)]);
    plot_field_robinson_replace(lon,lat,field,ctrs); %pcontinents; set(gca,'xlim',[0 357.5]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
    if max(lat)==86
        set(gca,'ylim',[-86 86])
    end
    colorbar off; originalSize = get(gca, 'Position'); hc = colorbar('eastoutside'); 
    if max(ctrs)==1
        hc = colorbar; set(hc,'ytick',-1:0.5:1); set(hc,'yticklabel',{'-1','-0.5','0','0.5','1�C'})
    elseif max(ctrs)==2
        hc = colorbar; set(hc,'ytick',-2:1:2); set(hc,'yticklabel',{'-2','-1','0','1','2�C'})
    elseif max(ctrs)==0.2
        set(hc,'ytick',-0.2:0.1:0.2)
    elseif max(ctrs)==0.6
        %caxis([-0.6 0.60001]); 
        set(hc,'ytick',-0.6:0.2:0.6)
    elseif max(ctrs)==0.8
        set(hc,'ytick',-0.8:0.4:0.8)
    end
    %set(gca,'xticklabel',''); set(gca,'yticklabel',''); 
    set(gca,'fontsize',16)
    set(gcf,'position',[0 1000 600 450]);
    originalSize(1) = originalSize(1)-0.04;
    set(gca, 'Position', originalSize);
    if nargin > 9
        title(title_text,'fontsize',14)
    end
    subplot(3,1,3)
else
    figure; pretty_figure(500,250,'Model Year','Standard Deviations',0:1200:6000,-10:2:10,16,0:100:500,-10:2:10);
end    

n = size(Xdisc,1);
ne = n/l;

dt = 1;

Xdisc_mean = zeros(l,1);

is_include = 1:ne;

t = (0:l-1)*dt+year1;

for i = is_include(end):-1:is_include(1)
    i1 = (i-1)*l+1;
    i2 = i*l;
    hold on;
    plot(t,Xdisc(i1:i2,ndisc),'color',[0.7 0.7 0.7]);
    Xdisc_mean = Xdisc_mean + Xdisc(i1:i2,ndisc);
end

i1 = (ne-1)*l+1;
i2 = ne*l;
plot((0:l-1)*dt+year1,Xdisc(i1:i2,ndisc),'linewidth',1.5);
%     
%     i1 = (imax-1)*l+1;
%     i2 = imax*l;
%     plot((0:l-1)*dt+year1,Xdisc(i1:i2,ndisc),'color',[0.7 0.2 0],'linewidth',1.5);
% end

Xdisc_mean = Xdisc_mean./length(is_include);
plot(t,Xdisc_mean,'k','linewidth',1.5);

originalSize = get(gca, 'Position');
originalSize(1) = originalSize(1)-0.04;
set(gca, 'Position', originalSize);

set(gca,'xlim',[floor(t(1)) ceil(t(end))]);
set(gca,'ylim',[floor(min(Xdisc(:,ndisc))) ceil(max(Xdisc(:,ndisc)))])
set(gca,'fontsize',14)
ylabel('Standard Deviations','fontsize',14);
%xlabel('Year','fontsize',14);