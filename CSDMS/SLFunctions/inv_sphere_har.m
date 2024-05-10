function [ F_yx ] = inv_sphere_har( vec_in, maxdeg, N, varargin)

% if not already generated, generate Legendre Polynomials
if nargin == 3
    [xGQ,~] = GaussQuad(N);
    x_GL = acos(xGQ)*180/pi - 90;
    colat_rad = (90-x_GL)*pi/180;
    x = cos(colat_rad);
    for l=0:N
        P_lm{l+1} = legendre(l,x,'norm');
    end
else
    P_lm = varargin{1};
end

% sum over degree and order
F_ym = zeros([1 2]* (N));
for l=0:maxdeg
    F_ym(:,1:l+1) = F_ym(:,1:l+1) + P_lm{l+1}(1:l+1,:).' .* repmat(get_coeffs(vec_in,l),N,1);
end

% Use inverse fft to sum over exponential
F_ym(:,N+2:end) = conj(F_ym(:,N:-1:2));
F_yx = ifft(F_ym,[],2);

% additional factors
F_yx = F_yx*(2*N)*sqrt(2);
% *(length of data) = *(2*maxdeg) needed for ifft in matlab
% *sqrt(2) as P_00 = 1/sqrt(2) for legendre(l,x,'norm')
end