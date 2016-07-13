function [ group ] = NextLevelStructure( CM ,clsIDX,thr)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 1
    thr = 333;
end
A = (CM+CM')/2;
Ncls = size(CM,1);
for i = 1:Ncls
    group{i} =clsIDX( i);
end
A = ones(Ncls)./(A+0.00001);
tmp = 1;
dis = zeros(Ncls*(Ncls-1)/2,1);
for i = 1: Ncls
    dis(tmp:tmp+Ncls-i-1) = A(i,i+1:end);
    tmp = tmp+ Ncls -i ;
end
Z = linkage(dis','average');
i = 1;
while(i<=size(Z,1))
    if Z(i,3)>thr
        break;
    end
    i = i+1;
end
N =Ncls + i-1;
for i = Ncls+1:N
    group{i} = [group{Z(i-Ncls,1)},group{Z(i-Ncls,2)}];
end
visit = zeros(Ncls,1);
idx = 1;i=0;
flag = 0;
reserve = [];
while(~isempty(find(visit == 0)))
    index = group{N-i};
    if visit(backTo(clsIDX,index(1))) == 0
        for j = 1:length(index)
            visit(backTo(clsIDX,index(j))) = 1;
        end
        reserve = [reserve,N-i];        
    end
    i = i+1;
end
group = group(reserve);
end
    
function[cls] = backTo(clsIDX,label)
for cls = 1:length(clsIDX)
    if(clsIDX(cls) == label)
        break;
    end
end
end
    


