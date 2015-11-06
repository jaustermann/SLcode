function T_lm = get_tlm(maxdeg)

a = 6371E3;
M_e = 5.9742E24;
const = 4*pi*a^3/M_e;

T_lm = [];
T = zeros(length(maxdeg));

for n = 0:maxdeg
    T(n+1) = const/(2*n + 1);
    T_add = repmat(T(n+1),1,n+1);
    T_lm = [T_lm T_add];
end