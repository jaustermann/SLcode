function spa_matrix = sph2spa_old(a_lm,maxdegree,lon,colat)

% Function to expand spherical harmonics coefficients that are given in the
% complex form into spatial matrix. lon and colat (in vector form) are the
% points on which the harmonics are evaluated.
% J. Austermann 2012


% Initialize spatial matrix
spa_matrix = zeros(length(colat),length(lon));

% loop over all degrees
for n = 0:maxdegree
    % get associated legendre polynomial in normalized form
    P_lm = legendre_me(n,cos(colat'*pi/180),'me');
    % get coefficients corresponding to degree n from coefficient vector
    a_l = get_coeffs(a_lm,n);
        
    % initialize m vector and evaluate exponential part of spherical
    % harmonic
    m_vec = 0:n;
    exp_func = exp(1i*m_vec'*lon*pi/180);
    
    % add / sum all orders from m = 0 to m = n to spatial matrix
    % m = 0
    Mat_sum = (P_lm(1,:)'*exp_func(1,:)) * a_l(1);
    spa_matrix = spa_matrix + Mat_sum;
    
    % m = 1:n
    for m = 1:n
        k = m+1;
        Mat_sum = 2 * real((P_lm(k,:)'*exp_func(k,:)) * a_l(k));
        spa_matrix = spa_matrix + Mat_sum;
    end
    
end
        
end