function [D] = build_dict(trainsplit,validsplit,clsIDX,fromD,fromvalid)
global option 
D = [];
if nargin == 3
    n = length(clsIDX);
    fromD = cell(n,1);
    fromvalid = cell(n,1);
end
for cls = 1:length(clsIDX)
    D = [D,fromD{clsIDX(cls)}];
    n = length(fromvalid{clsIDX(cls)});
    D = [D,validsplit{clsIDX(cls)}(:,fromvalid{clsIDX(cls)}(1:min(n,option.M/10)))];
    n = size(fromD{clsIDX(cls)},2)+min(n,option.M/10);
    r = option.M-n;
    idx = randperm(size(trainsplit{clsIDX(cls)},2),r);
    D = [D,trainsplit{clsIDX(cls)}(:,idx)];
end
end