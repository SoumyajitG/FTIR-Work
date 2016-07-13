function [correct,wrong]=evaluate(Y,label,D,W,clsIDX,P)
global option
load_options;
correct = 0;wrong = 0;
for i = 1:size(Y,2)
    idx=1;
    while idx<= length(D)
        s=[];cls=[];
        SC = OMP(D{idx},Y(:,i),option.T);
        
        for j = 1 : option.H + 1
            s(:,j) = W{idx,j}*P{idx,j}*SC;
        end
        cls(:,i) = vote1(s);
        c = sum(clsIDX{idx}*cls(:,i));
        j = idx+1;
        for j = idx+1:length(clsIDX)
            if ~isempty(find(clsIDX{j} == c))
                idx = j;
                break;
            end
        end
        if (j==length(clsIDX))
            if  (idx<j)
                break;
            end
        end
    end
    pred(i) = c;
    gt(i) = sum((1:10)*label(:,i));
    if pred(i) == gt(i)
        correct = correct + 1;
    else
        wrong = wrong+1;
    end
end
end

