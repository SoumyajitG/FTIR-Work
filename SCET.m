clear variables;
global option
load_options;
 [trainsplit,validsplit,testsplit] = load_data_mnist;
 clsIDX = 1:10;
%load('SCET.mat');
P = gen_P(length(clsIDX));
tmp = 0;
totalP = P;
CLSIDX{1} = clsIDX;
fromD = cell(10,1);fromvalid = cell(10,1);
while (size(totalP,1) >  tmp&&tmp<8)
    tmp = tmp+1;
    D = build_dict(trainsplit,validsplit,CLSIDX{tmp},fromD,fromvalid);
    [train,valid,trainlebal,validlebal,trainidx,valididx] = prepare_train_data_mnist(trainsplit,validsplit,CLSIDX{tmp});
    [group,W,DICTIONARY{tmp},fromvalid,fromD] = singleLevel(D,train,valid,trainlebal,validlebal,CLSIDX{tmp},totalP(tmp,:),valididx);
    CLASSIFIER(tmp,:) = W;
    for i = 1:length(group)
        if (length(group{i})>1)
            totalP = [totalP;gen_P(length(group{i}))];
            CLSIDX = [CLSIDX;group(i)];
        end
    end
end
[test,testlebal]= prepare_test_data_mnist( testsplit);
[correct,wrong] = evaluate(test, testlebal, DICTIONARY, CLASSIFIER, CLSIDX, totalP);



    
    
    
  



