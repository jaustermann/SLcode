% Script to read in and save the default topography from etopo

cd '/Users/jackyaustermann/Documents/MATLAB'
[Z, refvec] = etopo(15);
cd '/Users/jackyaustermann/Desktop/SLcode/SLFunctions'

lon_topo = linspace(-180,180,length(Z(1,:))+1);
lon_topo = lon_topo(1:end-1);
lon_topo = lon_topo+(lon_topo(2)-lon_topo(1))/2;

% change from -180 180 to 0 360
topo_orig = [Z(:,end/2:end) Z(:,1:end/2+1)];
lon_topo = [lon_topo(end/2:end) lon_topo(1:end/2+1)+360];

lat_topo = linspace(-90,90,length(Z(:,1))+1);
lat_topo = lat_topo(1:end-1);
lat_topo = lat_topo +(lat_topo(2)-lat_topo(1))/2;

save('topo_SL','lat_topo','lon_topo','topo_orig')