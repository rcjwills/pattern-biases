
load('tk_obs.mat')
tk_obs = -tk_obs(:,1:2,:,:);

ndisc = 1; [tk1_mean,tk1_anom] = plot_multivariate_model_bias_patterns(-tk,-FPs,ndisc,42,1979,LON_AXIS,LAT_AXIS,sst_scale,slp_scale,linspace(-0.4,0.4,25),linspace(-80,80,25));
ndisc = 2; [tk2_mean,tk2_anom] = plot_multivariate_model_bias_patterns(-tk,-FPs,ndisc,42,1979,LON_AXIS,LAT_AXIS,sst_scale,slp_scale,linspace(-0.4,0.4,25),linspace(-80,80,25));

for i = 1:size(tk1_anom,2)
    p = polyfit((1:42)',tk1_anom(:,i),1);
    tk1_trends(i) = p(1).*41;
    p = polyfit((1:42)',tk2_anom(:,i),1);
    tk2_trends(i) = p(1).*41;
end

for i = 1:3
    for j = 1:2
        p = polyfit((1:42)',tk_obs(:,1,i,j),1);
        tk1_obs_trends(i,j) = p(1)*41;
        p = polyfit((1:42)',tk_obs(:,2,i,j),1);
        tk2_obs_trends(i,j) = p(1)*41;
    end
end