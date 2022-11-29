
function RSL_spec = get_spectrum(RSL_lm, maxdeg)
RSL_pow = zeros(1,maxdeg+1);

ind = 1;
for j = 1:maxdeg+1
    for i = 1:j
        if i == 1
            RSL_pow(j) = RSL_pow(j) + abs(RSL_lm(ind))^2;
        else
            RSL_pow(j) = RSL_pow(j) + 2*abs(RSL_lm(ind))^2;
        end
        ind = ind+1;
    end
end

%l_deg = 0:maxdeg;
RSL_spec = RSL_pow * 4*pi; %./(2*l_deg+1); Don't need this for fully normalized coefficients
end

%figure
%plot(0:1024,abs(RSL_pow))