function [test,testlebal]= prepare_test_data_mnist( testsplit)
global option
test=[];testlebal=[];
for i = 1:10
    rd = randperm(size(testsplit{i},2),100);
        test = [test,testsplit{i}(:,rd)];
      
        l = zeros(10,1);l(i) = 1;
        testlebal = [testlebal,repmat(l,1,size(testsplit{i}(:,rd),2))];
    
end