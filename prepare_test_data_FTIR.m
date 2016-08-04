function [test,testlebal]= prepare_test_data_FTIR( testsplit)
global option
test=[];testlebal=[];
for i = 1:16
        test = [test,testsplit{i}];
      
        l = zeros(16,1);l(i) = 1;
        testlebal = [testlebal,repmat(l,1,size(testsplit{i},2))];
    
end