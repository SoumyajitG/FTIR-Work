function [cm,SC,fromvalid] = classification(D,Y,W,P,clsIDX,L,valididx)
global option
Nsample = size(Y,2);
%SC = iter_SC(D,L,Y,W,P);
for i = 1:Nsample
   SC(:,i) = OMP(D,Y(:,i),option.T);
%SC(:,i) = OMP(D,Y(:,i),30);
for j = 1 : option.H+1
    s(:,j) = W{j}*P{j}*SC(:,i);
end
    cls(:,i) = vote1(s);
end
[cm,fromvalid] = computeCM(cls,L,clsIDX,valididx);
end

function [cm,fromvalid] = computeCM(cls,L,clsIDX,valididx)
global option
fromvalid = cell(10,1);
Ncls = size(L,1);
Nsample = size(L,2);
cm = zeros(Ncls);

for i = 1: Nsample
    x = (1:Ncls)*L(:,i);
    y = (1:Ncls)*cls(:,i);
    clsx = clsIDX*L(:,i);
    clsy = clsIDX*cls(:,i);
%     if x~= y
%         tmp = mod(i,option.validNumPerCls);
%         if tmp == 0
%             tmp = tmp+option.validNumPerCls;
%         end
%         fromvalid{clsx} = [fromvalid{clsx};valididx{clsx}(tmp)];
%     end
    cm(y,x) = cm(y,x) + 1;
end
for i = 1:Ncls
    cm(:,i) = cm(:,i)/sum(cm(:,i));
end
end

