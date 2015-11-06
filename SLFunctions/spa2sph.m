function a_lm = spa2sph(func_gauss,maxdegree,lon,colat,varargin)

% Function to decompose spatial values on a sphere (given lon and lat in 
% vector form) into spherical harmonic coefficients in the complex form. 
% It follows the concept of fast spherical harmonic decomposition of 
% 'Numerical recipes' of Press et al.
% J. Austermann 2012

% Obtain weights and points for gaussian quadrature. 
% Takes as many points as they are given in lat direction
N = length(colat);
[x,w] = GaussQuad(N);
% The argument of Legendre polynomials (cos(x)) are quadrature points,
% hence one has to calculate the respective latitude points
x_gauss = acos(x)*180/pi;

% % check if you already are on a gauss legendre grid
% if sum(abs(x_gauss - colat)) < 0.1
%     func_gauss = func;
% else
%     % Interpolate to gaussian quadrature points in lat direction
%     [lon_out, lat_out] = meshgrid(lon,x_gauss);
%     [lon_in, lat_in] = meshgrid(lon,colat);
%     func_gauss = interp2(lon_in, lat_in, func, lon_out, lat_out);
% end
    
% Sum over longitude (can also be done using fft)

% g_m = zeros(N,maxdegree+1);
% for m = 0:maxdegree
%     exp_vec = exp(-1i*m*lon*pi/180);
%     exp_mat = repmat(exp_vec,N,1);
%     vec = func_gauss.*exp_mat;
%     g_m(:,m+1) = sum(vec,2);
% end

m = 0:maxdegree;
g_m = func_gauss * exp(-1i*lon'*m*pi/180);

% scale matrix by longitude increment (only works for equally spaced
% grids!)
del_lon = lon(2) - lon(1);
g_m = g_m*del_lon*pi/180;


% Do Gaussian quadrature in latitude direction. Legendre_me is a
% normalization of legendre polymonials according to JXM.

% check whether the Legendre Polynomials have been precomputed
if nargin == 4

    % if not compute them here
    a_lm = [];
    for n = 0:maxdegree
        % Get legendre polynomial
        P_lm(1:n+1,1:N) = legendre_me(n,cos(x_gauss*pi/180),'me');

        % Do quadrature
        a_lm = [a_lm  1/(4*pi) * w*(P_lm'.*g_m(:,1:n+1))];
    end
    
    % if the have use them for the quadrature
else
    P_lm_spa2sph = varargin{1};
    
    a_lm = [];
    ind = 1;
    for n = 0:maxdegree
        % Get legendre polynomial
        for m = 0:n
            P_lm(m+1,1:N) = P_lm_spa2sph(ind:ind+N-1); 
            ind = ind + N;
        end
        % Do quadrature
        a_lm = [a_lm  1/(4*pi) * w*(P_lm'.*g_m(:,1:n+1))];
    end

end



end
