
function a_n = get_coeffs(a_lm,n)

if n == 0;
    a_n = a_lm(1);
else
    
    gauss_sum = n*(n+1)/2;
    % position of current order in a_lm vector
    vec = gauss_sum+1:gauss_sum+1+n;
    a_n = a_lm(vec);

end
    
end