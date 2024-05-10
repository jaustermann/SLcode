function [P_lm_spa2sph, P_lm_sph2spa] = get_Legendre(lat_in,maxdeg)
% Script to precompute Legendre polynomials so that they don't have to be
% computed each time
% J. Austermann 2015

% get latitude and colatitude
% lat = lat_in;
colat = 90 - lat_in;
%colat_inv = lat_in + 90;

% make a legendre grid
% N = length(lat);
% [x,w] = GaussQuad(N);
% The argument of Legendre polynomials (cos(x)) are quadrature points,
% hence one has to calculate the respective latitude points
%x_gauss = colat_inv; %acos(x)*180/pi;

% initialize matrices
P_lm_sph2spa = zeros(length(lat_in),maxdeg+1,maxdeg+1);
P_lm_spa2sph = zeros(maxdeg+1,length(lat_in),maxdeg+1); 

for n = 0:maxdeg
    % compute Legendre polynomials
    % spa2sph - needs to be on GL grid
    P_lm_spa2sph_n = sqrt(2)*legendre_me(n,cos(colat*pi/180),'norm');
    % sph2spa - can be on any grid
    P_lm_sph2spa_n = sqrt(2)*legendre_me(n,cos(colat'*pi/180),'norm');

    % write them int a vector
    for m = 0:n
        P_lm_sph2spa(1:length(lat_in),m+1,n+1) = P_lm_sph2spa_n(m+1,:)';
        P_lm_spa2sph(m+1,1:length(lat_in),n+1) = P_lm_spa2sph_n(m+1,:);
    end
    
end

end