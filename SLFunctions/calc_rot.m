function La_lm = calc_rot(L_lm,k_L,k_T)

% constants:
M = 5.976E24;
omega = 7.292e-5;
% g = 9.80665 ;
a = 6371e3;
k_hydro = 0.934;
% k_hydro = 0.934055;
CminA = 2.6E35;

% extract degree 2:
L20 = L_lm(4);
L21 = L_lm(5); 
L22 = L_lm(6);

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

% equations after Jerry's notes

I13 = sqrt(10/3) * a * M * real(L21);
I23 = - sqrt(10/3) * a * M * imag(L21);

% equation from Mitrovica and Wahr 2005 page 2&3

II = I13+1i*I23;

m0 = II/CminA *(1+k_L)/(1-(k_T/k_hydro));

m1 = real(m0);
m2 = imag(m0);
m3 = 0;

m = [m1, m2, m3];

% calculate the perturbation to the rotational potential from Mitrovica et
% al. 2001 page 3

La00 = a^2*omega^2/3*(m*m' + 2*m3);
La20 = a^2*omega^2/(6*sqrt(5))*(m1^2 + m2^2 - 2*m3^2 - 4*m3);
La21 = a^2*omega^2/sqrt(30)*(m1*(1+m3) - 1i*m2*(1+m3));
La22 = a^2*omega^2/(sqrt(5)*sqrt(24))*((m2^2-m1^2)+1i*2*m1*m2);
% La2m1 = -1 * conj(La21);
% La2m2 = 1 * conj(La22);

La_lm = zeros(size(L_lm));

La_lm(1:6) = [La00 0 0 La20 La21 La22];

end