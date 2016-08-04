function [train,valid,trainlebal,validlebal,trainidx,valididx]= prepare_train_data_FTIR(trainsplit,validsplit,cls)
global option
tmp=1;
train = [];trainlebal = [];valid = [];validlebal = [];
for i = 1:16
    if ~isempty(find(cls == i))
        trainidx{i} = randperm(size(trainsplit{i},2),option.trainNumPerCls);
        valididx{i} = randperm(size(validsplit{i},2),option.validNumPerCls);
        train = [train,trainsplit{i}(:,trainidx{i})];
        valid = [valid,trainsplit{i}(:,valididx{i})];
        l = zeros(length(cls),1);l(tmp) = 1;
        trainlebal = [trainlebal,repmat(l,1,option.trainNumPerCls)];
        validlebal= [validlebal,repmat(l,1,option.validNumPerCls)];
        tmp =tmp+1;
    end
end
end