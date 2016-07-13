function [group,W,D,fromvalid,fromD] = singleLevel(D,train,valid,trainlebal,validlebal,clsIDX,P,valididx)
global option
thr = 88;
for Ns = 1:size(train,2)
    SC(:,Ns) = OMP(D,train(:,Ns),option.T);
end

for iter=1:5%option.iter
    W = iter_W(trainlebal,P,SC);
    SC = iter_SC(D,trainlebal,train,W,P);
    D = iter_dict(train,D,SC);
%     for Ns = 1:size(Y,2)
%         SC(:,Ns) = OMP(D,Y(:,Ns),option.T);
%     end
%    [cm,SC,misIDX,~,misINFO]=classification(D,Y,W,P,clsIDX,L);
   disp('bingo');
end
%M = SC;
%W = Ltmp*M'*inv(M*M'+option.gamma/option.beta*eye(size(M,1)));
[cm,SC,fromvalid]=classification(D,valid,W,P,clsIDX,validlebal,valididx);
group = NextLevelStructure( cm ,clsIDX,thr);
fromD = cell(10,1);
tmp = 1;
for i = 1:10
    if ~isempty(find(clsIDX == i))
        sctmp = SC((tmp-1)*option.M+1:tmp*option.M,:);
        this = mean(sctmp(:,(tmp-1)*option.validNumPerCls+1:tmp*option.validNumPerCls),2);
        all = mean(sctmp,2);
        score = this./all;
        [~,rank] = sort(score,'descend');
        index = (tmp-1)*option.M + (rank(1:option.M/2));
        fromD{clsIDX(tmp)} = D(:,index);
        tmp = tmp + 1;
    end
end
end


