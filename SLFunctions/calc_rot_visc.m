function [La_lm, sdelI, sdelm] = calc_rot_visc(L_lm,k_L,k_T,t_it,beta_konly_l, ...
    beta_konly_tide, sdelI, sdelm)

% constants:
G = 6.67408E-11;
omega = 7.292E-5;
a = 6371000;
CminA = 2.6347269E35;

k_f = (3*G*CminA)/(a^5*omega^2);

% extract degree 2:
L00 = L_lm(1);
L20 = L_lm(4);
L21 = L_lm(5); 

% WHY ONLY CONSIDERING K NOT ALSO 1+K-H


% LOOK WHERE THESE ARE FROM -> notes
I(1) = sqrt(32/15)*pi*a^4*real(L21);
I(2) = - sqrt(32/15)*pi*a^4*imag(L21);
I(3) = (8/3)*pi*a^4*(L00 - L20/sqrt(5));

% calculate viscous contribution
if t_it == 2
    V_lm = [0,0,0];
    V_lm_T = [0,0,0];
else 
    V_lm = beta_konly_l{t_it-1} * sdelI(1:t_it-2,:);
    V_lm_T = beta_konly_tide{t_it-1} * sdelm(1:t_it-2,:);
end

m = 1/(1-k_T/k_f)*(1/CminA * ((1+k_L)*I + V_lm) + V_lm_T/k_f);
m(3) = 0;

sdelI(t_it-1,:) = I - sum(sdelI(1:t_it-2,:),1);
sdelm(t_it-1,:) = m - sum(sdelm(1:t_it-2,:),1);

% calculate the perturbation to the rotational potential from Milne 1998

La00 = a^2*omega^2/3*(m*m' + 2*m(3));
La20 = a^2*omega^2/(6*sqrt(5))*(m(1)^2 + m(2)^2 - 2*m(3)^2 - 4*m(3));
La21 = a^2*omega^2/sqrt(30)*(m(1)*(1+m(3)) - 1i*m(2)*(1+m(3)));
La22 = a^2*omega^2/(sqrt(5)*sqrt(24))*((m(2)^2-m(1)^2)+1i*2*m(1)*m(2));

La_lm = zeros(size(L_lm));

La_lm(1:6) = [La00 0 0 La20 La21 La22];

end