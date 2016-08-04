function [group ] = NextLevelStructureSoft( CM ,clsIDX,thr)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 1
    thr = 0.1;
end
A = (CM+CM')/2;
Ncls = size(CM,1);
for i = 1:Ncls
    group{i} =clsIDX( i);
    for j = 1:Ncls
        if (i~=j&&(CM(i,j)>thr || CM(j,i)>thr))
            group{i} = [group{i},clsIDX(j)];
        end
    end
end
idx = [];
for i = 1:Ncls
    for j = 1:i-1
        if (sum(ismember(group{i},group{j}) )== length(group{j})) 
            idx = [idx;j];
        elseif (sum(ismember(group{j},group{i}) )== length(group{i}))
            idx = [idx;i];
        end
    end
end
group(idx) = [];
end
   

