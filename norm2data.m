function [d] = norm2data(data)
d = data./repmat(sqrt(sum(data.*data,1)),size(data,1),1);
end