function spa_matrix = sph2spa(a_lm,maxdegree,lon,colat,varargin)

% Function to expand spherical harmonics coefficients that are given in the
% complex form into spatial matrix. lon and colat (in vector form) are the
% points on which the harmonics are evaluated.
% J. Austermann 2012
% updated 2015


% Initialize spatial matrix
spa_matrix = zeros(length(colat),length(lon));

% initialize m vector and evaluate exponential part of spherical
% harmonic
% Need to use '-' to get the longitude right. 
m_vec = 0:maxdegree;
exp_func_all = exp(-1i*m_vec'*lon*pi/180);


% check whether the Legendre Polynomials have been precomputed
if nargin == 4
    
    % if not compute them here

    % loop over all degrees
    for n = 0:maxdegree
        % get associated legendre polynomial in normalized form
        P_lm = legendre_me(n,cos(colat*pi/180),'me');

        % get coefficients corresponding to degree n from coefficient vector
        a_l = get_coeffs(a_lm,n);
        a_l_rep = repmat(a_l',1,length(lon));

        % pick out right exp_func
        exp_func = exp_func_all(1:n+1,:);

        % add / sum all orders from m = 0 to m = n to spatial matrix
        % m = 0
        % Mat_sum = (P_lm(1,:)'*exp_func(1,:)) * a_l(1);
        spa_matrix = spa_matrix + (P_lm(1,:)'*exp_func(1,:)) * a_l(1);

        % m = 1:n
        spa_matrix = spa_matrix + 2 * real(P_lm(2:end,:)' * (exp_func(2:end,:) .* a_l_rep(2:end,:)));

    end
    
    
    % if the have use them for the quadrature
else
    P_lm_sph2spa = varargin{1};
    
    % loop over all degrees
    ind = 1;
    for n = 0:maxdegree
        % get associated legendre polynomial in normalized form
        for m = 0:n
            P_lm(m+1,1:length(colat)) = P_lm_sph2spa(ind:ind+length(colat)-1); 
            ind = ind + length(colat);
        end

        % get coefficients corresponding to degree n from coefficient vector
        a_l = get_coeffs(a_lm,n);
        a_l_rep = repmat(a_l',1,length(lon));

        % pick out right exp_func
        exp_func = exp_func_all(1:n+1,:);

        % add / sum all orders from m = 0 to m = n to spatial matrix
        % m = 0
        % Mat_sum = (P_lm(1,:)'*exp_func(1,:)) * a_l(1);
        spa_matrix = spa_matrix + (P_lm(1,:)'*exp_func(1,:)) * a_l(1);

        % m = 1:n
        spa_matrix = spa_matrix + 2 * real(P_lm(2:end,:)' * (exp_func(2:end,:) .* a_l_rep(2:end,:)));

    end
    
end

end