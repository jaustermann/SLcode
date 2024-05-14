
% This script plots how much sea level is changing from the first to the
% second timestep. This can be used to plot a sea level fingerprint. Here
% we normalize the sea level change so that the average sea level change
% over the ocean basins is 1. This allows us to better compare fingerprints
% from different melt scenarios. 

%% Plot results: fingerprint

% For fingerprint, look at change in RSL
del_SL = RSL(:,:,2) - RSL(:,:,1);

% calculate the scaling to normalize RSL (it's normalized to be
% one on average, when averaged over the final ocean basin).
FP_normalization_factor = oc0_lm(1)/sphere_har(del_SL.*oc_0,0,N,P_lm); 

% What to plot: (e.g. RSL, ice, topo)
plot_timeslice = squeeze(del_SL*FP_normalization_factor); 

figure
hold on
pcolor(lon_out,lat_out,plot_timeslice)

% plot today's coastlines
lon_shift = [lon_topo(end/2+1:end-1)-360 lon_topo(2:end/2)];
topo_shift = [topo_ice(:,end/2+1:end-1) topo_ice(:,2:end/2)];
contour(lon_shift,lat_topo,topo_shift,[0 0],'k')

shading flat
colorbar 

title('Normalized sea level change')