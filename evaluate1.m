function [correct,wrong]=evaluate1(Y,label,D,W,clsIDX,P,child)
global option
load_options;
correct = 0;wrong = 0;
for i = 1:size(Y,2)
    disp(i);
    cur = 1;list = 1;
    tmp = 1;class = [];
    while(tmp<=length(list))
        s=[];cls=[];
        SC = OMP(D{list(tmp)},Y(:,i),option.T);
        for j = 1 : option.H + 1
            s(:,j) = W{tmp,j}*P{tmp,j}*SC;
        end
        cls(:,i) = vote1(s);
        c = sum(clsIDX{tmp}*cls(:,i));
        flag = 0;
        if length(child) > 1 
            for j = child{cur}(1):child{cur}(end)
                if ismember(clsIDX{j},c)
                    list = [list;j];
                    flag = 1;
                end
            end
        end
        if flag ==0
           class = [class;c];
        end
        tmp = tmp+1;
    end
    pred(i) = argmax(class);
     gt(i) = sum((1:size(label,1))*label(:,i));
    if pred(i) == gt(i)
        correct = correct + 1;
    else
        wrong = wrong+1;
    end
            
                
end
end

