% Script to read in and save the default topography from etopo

cd '/Users/jackyaustermann/Documents/MATLAB'
[Z, refvec] = etopo(15);
cd '/Users/jackyaustermann/Desktop/SLcode/SLFunctions'

% change from -180 180 to 0 360
Z = [Z(:,end/2-1:end) Z(:,1:end/2)];
topo_orig = Z;

lon_topo = linspace(0,360,length(Z(1,:)));
lon_topo = lon_topo(1:end-1);
lon_topo = lon_topo+(lon_topo(2)-lon_topo(1))/2;
lon_topo = [lon_topo(end)-360 lon_topo];

lat_topo = linspace(-90,90,length(Z(:,1))+1);
lat_topo = lat_topo(1:end-1);
lat_topo = lat_topo +(lat_topo(1)/lat_topo(2))/2;

save('topo_SL','lat_topo','lon_topo','topo_orig')