function La_lm = calc_rot(L_lm,k_fl,k_fl_tide,length_of_day)

% constants:
% G = 6.67408E-11;
% M = 5.972E24;

if nargin == 4
    omega = 2*pi/length_of_day/(60*60);
else
    omega = 7.292E-5;
end

% g = 9.81;
a = 6371000;
%k_hydro = 0.934;
%of 100km and mantle convection support
% k_hydro = 0.934055;
CminA = 2.6347269E35;
%CminA = 2.6347269E35;

%k_hydro = (3*G*CminA)/(a^5*omega^2);
k_hydro = 0.942; % this is for a rotating earth with lithospheric thickness 

% extract degree 2:
L20 = L_lm(4);
L21 = L_lm(5); 
L22 = L_lm(6);
k_L = k_fl(2);
k_T = k_fl_tide(2);

% calculate ineartia tensor after Matsuyama et al. 2006, page 6

% I11 = 4*pi*a^4*(1+k_L)*(1/(3*sqrt(5))*L20 - sqrt(2/15)*real(L22)) - ...
%     M*omega^2*a^3/(9*g) * (k_hydro - k_T);
% 
% I22 = 4*pi*a^4*(1+k_L)*(1/(3*sqrt(5))*L20 + sqrt(2/15)*real(L22)) - ...
%     M*omega^2*a^3/(9*g) * (k_hydro - k_T);
% 
% I33 = -8*pi*a^4/(3*sqrt(5))*(1+k_L)*L20 + 2*M*omega^2*a^3/(9*g) * ...
%     (k_hydro - k_T);
% 
% I12 = 8*pi*a^4/sqrt(30)*(1+k_L)*imag(L22);
% 
% I13 = 8*pi*a^4/sqrt(30)*(1+k_L)*real(L21);
% 
% I23 = -8*pi*a^4/sqrt(30)*(1+k_L)*imag(L21);

I13 = sqrt(32/15)*pi*a^4*real(L21);
I23 = - sqrt(32/15)*pi*a^4*imag(L21);

II = I13+1i*I23;

% equation from Mitrovica and Wahr 2005

m0 = II/CminA *(1+k_L)/(1-(k_T/k_hydro));

m1 = real(m0);
m2 = imag(m0);
m3 = 0;

m = [m1, m2, m3];

% calculate the perturbation to the rotational potential from Milne 1998

La00 = a^2*omega^2/3*(m*m' + 2*m3);
La20 = a^2*omega^2/(6*sqrt(5))*(m1^2 + m2^2 - 2*m3^2 - 4*m3);
La21 = a^2*omega^2/sqrt(30)*(m1*(1+m3) - 1i*m2*(1+m3));
La22 = a^2*omega^2/(sqrt(5)*sqrt(24))*((m2^2-m1^2)+1i*2*m1*m2);
%La2m1 = -1 * conj(La21);
%La2m2 = 1 * conj(La22);

La_lm = zeros(size(L_lm));

La_lm(1:6) = [La00 0 0 La20 La21 La22];

end