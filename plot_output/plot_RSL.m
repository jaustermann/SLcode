
% This script plots relative sea level (RSL) change, i.e. sea level at a given
% time relative to today. The first part plots spatially varying RSL at a 
% specified time. The second part plots RSL at a specified location through 
% time. 

%% Plot results: time slice RSL

% Set time (in ka)
fig_time = 5; 
ind = find(ice_time_new==fig_time);

% What to plot: (e.g. RSL, ice, topo)
plot_timeslice = squeeze(RSL(:,:,ind)); 

figure
hold on
pcolor(lon_out,lat_out,plot_timeslice)

% plot coastlines at the given time
contour(lon_out,lat_out,topo(:,:,ind),[0 0],'w')

shading flat
colorbar 

title(['Relative sea level at ' num2str(fig_time) 'ka'])

%% Plot results: time series 

% Set latitude (-90 to 90) and longitude (-180 to 180)
lon_pt = 151.2;
lat_pt = -33.9;

clear RSL_pt
for i = 1:length(ice_time_new)
    RSL_pt(i) = interp2(lon_out,lat_out,squeeze(RSL(:,:,i)),lon_pt,lat_pt);
    %ice_thk_pt(i) = interp2(lon_out,lat_out,squeeze(ice(:,:,i)),lon_pt,lat_pt);
end

% plot
figure
hold on
plot(ice_time_new,RSL_pt','b')
plot(ice_time_new,GMSL_fo','r') % ice volume equivalent sea level spreading ice over a fixed (present-day) ocean area
plot(ice_time_new,GMSL_tvo','g') % ice volume equivalent sea level spreading ice over a time-varying ocean area
% plot(ice_time_new,ice_thk_pt') % plot thickness of ice sheet at given location
legend('Relative Sea Level','Global mean sea level (fixed ocean area)', 'Global Mean Sea Level (time varying ocean area)')
ylabel('Sea Level (m)')
xlabel('Time (ka)')

title('Sea level at a specified location')
