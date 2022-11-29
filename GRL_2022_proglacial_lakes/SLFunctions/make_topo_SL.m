% Script to read in and save the default topography from etopo

cd '~/Sync/MATLAB/etopo_ice_surf'
[Z, refvec] = etopo(15);

% change from -180 180 to 0 360
topo_ice = [Z(:,end/2:end) Z(:,1:end/2+1)];

cd '~/Sync/MATLAB/etopo_bed'
[Z, refvec] = etopo(5);
% change from -180 180 to 0 360
topo_bed = [Z(:,end/2:end) Z(:,1:end/2+1)];

% cell centered
lon_topo = linspace(-180,180,length(Z(1,:))+1);
lon_topo = lon_topo(1:end-1);
lon_topo = lon_topo+(lon_topo(2)-lon_topo(1))/2; % cell centered

% change from -180 180 to 0 360
lon_topo = [lon_topo(end/2:end) lon_topo(1:end/2+1)+360];

lat_topo = linspace(-90,90,length(Z(:,1))+1);
lat_topo = lat_topo(1:end-1);
lat_topo = lat_topo +(lat_topo(2)-lat_topo(1))/2; % cell centered


cd '/Users/jackyaustermann/Desktop/SLcode/SLFunctions'

save('topo_SL','lat_topo','lon_topo','topo_bed','topo_ice')

% write topo into ascii file
