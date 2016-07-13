function[P] = gen_P(ncls)
global option
P{1} = eye(option.M*ncls)
for i = 2:option.H+1
    p = zeros(option.numSelPerCls*ncls,option.M*ncls);
    tmp = 1;
    for j = 1:ncls
        idx = randperm(option.M,option.numSelPerCls);
        for k = 1:option.numSelPerCls
            p(tmp,(j-1)*option.M+idx(k)) = 1;
            tmp = tmp +1;
        end
    end
    P{i} = p;
end