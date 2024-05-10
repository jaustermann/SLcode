function func_GL = interpGL(func,lon,lat)

% Obtain weights and points for gaussian quadrature. 
% Takes as many points as they are given in lat direction
N = length(lat);
[x,w] = GaussQuad(N);
% The argument of Legendre polynomials (cos(x)) are quadrature points,
% hence one has to calculate the respective latitude points
x_gauss = acos(x)*180/pi;

% Interpolate to gaussian quadrature points in lat direction
[lon_out lat_out] = meshgrid(lon,x_gauss);
[lon_in lat_in] = meshgrid(lon,lat);
func_GL = interp2(lon_in, lat_in, func, lon_out, lat_out);

end