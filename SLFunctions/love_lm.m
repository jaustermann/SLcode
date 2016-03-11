function h_lm = love_lm(num, maxdegree)

% no dc shift also for gravity love number
h = [0 num'];

h_lm = [];

for n = 0:maxdegree
    
    h_add = repmat(h(n+1),1,n+1);
    h_lm = [h_lm h_add];
    
end