function [cls] = vote(score)
V = zeros(size(score));
for i = 1:size(score,2)
    idx = find(score(:,i) == max(score(:,i)));
    if length(idx)==0
        idx = 1;
    end
    V(idx(1),i) = 1;
end
V = sum(V,2);
id = find(V==max(V));
cls = zeros(size(score,1),1);
cls(id(1)) = 1;
end