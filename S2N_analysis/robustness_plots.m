models = {'access','canesm2','canesm5','cesm1','cesm2','cnrm_cm6','csiro_mk36', ...
'ec-earth3','gfdl_cm3','gfdl_esm2m','giss_e21g','ipsl_cm6a','miroc6', ...
'miroc_esm2l','mpi','norcpm1'};

load('16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N.mat', 's','tk','FPs')
tk1_emean = plot_multivariate_model_bias_patterns(-tk,-FPs,1,42,1979);
tk2_emean = plot_multivariate_model_bias_patterns(-tk,-FPs,2,42,1979);
p = polyfit(1:42,tk1_emean',1);
trend1 = 41*p(1);
p = polyfit(1:42,tk2_emean',1);
trend2 = 41*p(1);

for i = 1:16
    load(['16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N_robustness_',models{i},'_single_model.mat'],'s_single','tk_single')
    load(['16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N_robustness_',models{i},'_other_models.mat'],'s_other','tk_other')
    s_other_all(:,(i-1)*5+1:i*5) = s_other;
    s_single_all(:,(i-1)*5+1:i*5) = s_single;
    for j = 1:5
        tk1_single_emean = plot_multivariate_model_bias_patterns(tk_single(:,:,j),FPs,1,42,1979);
        tk2_single_emean = plot_multivariate_model_bias_patterns(tk_single(:,:,j),FPs,2,42,1979);
        tk1_other_emean = plot_multivariate_model_bias_patterns(tk_other(:,:,j),FPs,1,42,1979);
        tk2_other_emean = plot_multivariate_model_bias_patterns(tk_other(:,:,j),FPs,2,42,1979);
        p = polyfit(1:42,tk1_single_emean',1);
        trend1_single(j) = abs(41*p(1));
        p = polyfit(1:42,tk2_single_emean',1);
        trend2_single(j) = abs(41*p(1));
        p = polyfit(1:42,tk1_other_emean',1);
        trend1_other(j) = abs(41*p(1));
        p = polyfit(1:42,tk2_other_emean',1);
        trend2_other(j) = abs(41*p(1));
    end
    trend1_other_all(:,(i-1)*5+1:i*5) = trend1_other;
    trend1_single_all(:,(i-1)*5+1:i*5) = trend1_single;
    trend2_other_all(:,(i-1)*5+1:i*5) = trend2_other;
    trend2_single_all(:,(i-1)*5+1:i*5) = trend2_single;
end

load('16model_multivariate_ERSST5_ERA5_1979-2020_Amean_80deg_5yr_N20_S2N_robustness_random_crossmodel_sampling.mat', 's_random')
tmp = s_random(1,:);
tmp(31:80) = nan;

%%

figure; violinplot([tmp; s_single_all(1,:); s_other_all(1,:)]','','Width',0.4);
hold on; plot(linspace(4,8.75,20),s(1:20),'ko','markersize',7,'linewidth',2)
hold on; plot([0 10],[0.0282 0.0282],'k--','linewidth',1)
set(gca,'xlim',[1.4 9])
pretty_figure(300,400,'Pattern #','Signal Fraction',[1 2 3,linspace(4,8.75,20)],'none',16,{'','','','1','','3','','5','','7','','9','','','12','','','15','','','','','20'},'none')

figure; violinplot([tmp./(1-tmp); s_single_all(1,:)./(1-s_single_all(1,:)); s_other_all(1,:)./(1-s_other_all(1,:))]','','Width',0.4);
hold on; plot(linspace(4,8.75,20),s(1:20)./(1-s(1:20)),'ko','markersize',7,'linewidth',2)
hold on; plot([0 10],[0.029 0.029],'k--','linewidth',1)
set(gca,'xlim',[1.4 9])
pretty_figure(300,400,'Pattern #','Signal-to-Noise Ratio',[1 2 3,linspace(4,8.75,20)],'none',16,{'','','','1','','3','','5','','','','','10','','','','','15','','','','','20'},'none')

figure; violinplot([trend1_single_all; trend1_single_all; trend1_other_all]','','Width',0.4);
hold on; plot(4,trend1,'ko','linewidth',1.5,'markersize',7,'linewidth',2)
set(gca,'xlim',[1.4 4.2])
pretty_figure(175,400,'Pattern #','Trend',[1 2 3,linspace(4,8.75,20)],'none',16,{'','','1','','','3','','5','','','','','10','','','','','15','','','','','20'},'none')

figure; violinplot([trend2_single_all; trend2_single_all; trend2_other_all]','','Width',0.4);
hold on; plot(4,trend2,'ko','linewidth',1.5,'markersize',7,'linewidth',2)
set(gca,'xlim',[1.4 4.2])
pretty_figure(175,400,'Pattern #','Trend',[1 2 3,linspace(4,8.75,20)],'none',16,{'','','2','','','3','','5','','','','','10','','','','','15','','','','','20'},'none')

figure; violinplot([sqrt(trend1_single_all.^2+trend2_single_all.^2); sqrt(trend1_single_all.^2+trend2_single_all.^2); sqrt(trend1_other_all.^2+trend2_other_all.^2)]','','Width',0.4);
hold on; plot(4,sqrt(trend1.^2+trend2.^2),'ko','linewidth',1.5,'markersize',7,'linewidth',2)
set(gca,'xlim',[1.4 4.2])
pretty_figure(150,400,'','Combined Patterns 1 and 2 Trend Amplitude',[0 1],'none',16)
ylabel('Combined Patterns 1 and 2 Trend Amplitude')

lines6 = lines(6);
figure; h1 = plot(trend1_single_all,trend2_single_all,'o','color',lines6(6,:),'markerfacecolor',lines6(6,:),'markersize',5);
alpha(h1,.5)
hold on; h2 = plot(trend1_other_all,trend2_other_all,'o','color',lines6(4,:),'markerfacecolor',lines6(4,:),'markersize',5);
alpha(h2,.5)
hold on; plot(trend1,trend2,'ko','markersize',7,'linewidth',2)
set(gca,'xlim',[0 3.5])
pretty_figure(400,400,'Pattern-1 Trend','Pattern-2 Trend',0:0.5:4,0:0.5:4,16)

% figure; violinplot([trend1_single_all+trend2_single_all; trend1_single_all+trend2_single_all; trend1_other_all+trend2_other_all]','','Width',0.4);
% hold on; plot(4,trend1+trend2,'ko','linewidth',1.5)
% set(gca,'xlim',[1.5 4.2])
% pretty_figure(150,500,'Pattern #','Signal-to-Noise Ratio',[1 2 3,linspace(4,8.75,20)],'none',16,{'','','','2','','3','','5','','','','','10','','','','','15','','','','','20'},'none')