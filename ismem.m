function [bo] = ismem(clsidx,g)
n = length(clsidx);
bo = 0;
for i = 1:n
    nn = length(clsidx{i});
    p = sum(ismember(g,clsidx{i}));
    if p == nn
        bo = 1;
    end
end
        