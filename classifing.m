function [cls] = classifing(SC)
global option
for i = 1:size(SC,2)
    sc =reshape(SC(:,i),option.M,[]);
    sc = sum(sc,1);
    idx = find(sc==max(sc));
    cls(i) = idx(1);
end
end
    