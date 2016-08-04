clear variables;
global option
load_options;
% [trainsplit,validsplit,testsplit] = load_data_mnist;
%[trainsplit,validsplit,testsplit] = load_data_scene;
[trainsplit,validsplit,testsplit] = load_data_FTIR;
ncls = 16;
 clsIDX = 1:ncls;
%load('SCET.mat');
P = gen_P(length(clsIDX));
tmp = 0;
totalP = P;
CLSIDX{1} = clsIDX;
fromD = cell(ncls,1);fromvalid = cell(ncls,1);
for i = 1:length(trainsplit)
        r= randperm(size(trainsplit{i},2),option.M);
        fromD{i} = trainsplit{i}(:,r);
end

    D{1} = build_dict(trainsplit,validsplit,clsIDX,fromD,fromvalid);
%    D{1}  = rand(784,600);
%    D{1} = norm2data(D{1});
while (size(totalP,1) >  tmp&&tmp<12)
    
    tmp = tmp+1;
%    D = build_dict(trainsplit,validsplit,CLSIDX{tmp},fromD,fromvalid);
    [train,valid,trainlebal,validlebal,trainidx,valididx] = prepare_train_data_FTIR(trainsplit,validsplit,CLSIDX{tmp});
    [group,W,DICTIONARY{tmp},fromvalid,fromD] = singleLevel(D{tmp},train,valid,trainlebal,validlebal,CLSIDX{tmp},totalP(tmp,:),valididx);
    CLASSIFIER(tmp,:) = W;n = length(D);
     child{tmp} = [];
    for i = 1:length(group)
        if (length(group{i})>1&&(ismem(CLSIDX,group{i})==0 ))
            totalP = [totalP;gen_P(length(clsIDX))];
          % totalP = [totalP;gen_P(10)];
            CLSIDX = [CLSIDX;group(i)];
            D =[D, DICTIONARY(tmp)];
            n = length(D);
            child{tmp} = [child{tmp} ;n];
        end
    end
end
 [test,testlebal]= prepare_test_data_FTIR( testsplit);
[correct,wrong] = evaluate1(test, testlebal, DICTIONARY, CLASSIFIER, CLSIDX, totalP,child);



    
    
    
  



