function [correct,wrong]=evaluate(Y,label,D,W,clsIDX,P,child)
global option
load_options;
correct = 0;wrong = 0;
for i = 1:size(Y,2)
    disp(i);
    idx=1;
    while idx<= length(D)
        s=[];cls=[];
        SC = OMP(D{idx},Y(:,i),option.T);
        
        for j = 1 : option.H + 1
            s(:,j) = W{idx,j}*P{idx,j}*SC;
        end
        cls(:,i) = vote1(s);
        c = sum(clsIDX{idx}*cls(:,i));
        for j = child{idx}(1): child{idx}(end)
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

