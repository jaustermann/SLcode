function [ a_lm ] = sphere_har( F_yx, maxdeg, N, varargin)

% if not already generated, generate Legendre Polynomials
[x,w] = GaussQuad(N);
if nargin == 3
    for l=0:N
        P_lm{l+1} = legendre(l,x,'norm');
    end
else
    P_lm = varargin{1};
end

% Use fft to sum over exponential
F_ym = fft(F_yx,[],2);

% do Gauss Legendre quadrature
ind_a = 1;
a_lm = zeros(1,((maxdeg+1)*(maxdeg+2))/2);
for l=0:maxdeg
    a_lm(ind_a:ind_a+l) = sum(P_lm{l+1}(1:l+1,:).'.*repmat(w',1,length(1:l+1)).*F_ym(:,1:l+1),1);
    ind_a = ind_a+l+1;
end

% additional factors
a_lm = a_lm/(2*N)/sqrt(2);
% *(length of data) = *(2*maxdeg) needed for ifft in matlab
% *sqrt(2) as P_00 = 1/sqrt(2) for legendre(l,x,'norm')
end