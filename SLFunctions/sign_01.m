function out = sign_01(in)

out = -0.5*sign(in)+0.5;
out = 0.5*sign(out-0.6)+0.5;

end