function [ train,label,test,gt ,trainsplit,trainlen,D] = data_handler(  )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global option
gt=[];
train = [];label = [];D=[];test=[];
for i = 0:9
    data = load(['mnist/train',num2str(i),'.mat']);
    testdata =  load(['mnist/test',num2str(i),'.mat']);
    testdata.D = testdata.D(1:50,:);
    data.Dori = data.D;
    data.D = data.D(1:100,:);
    trainsplit{i+1} = data.D';
    trainlen(i+1) = size(data.D,1);
    idx = randperm(size(data.Dori,1),option.M);
    D = [D,data.Dori(idx,:)'];
    l = zeros(10,1);
    l(i+1) = 1;
    label = [label,repmat(l,1,size(data.D,1))];
    test = [test,testdata.D'];
    gt = [gt,repmat(l,1,size(testdata.D,1))];
end
for i = 1:10
    train = [train,trainsplit{i}];
end
train = train/255;
test = test/255;
D = D./sqrt(repmat(sum(D.*D,1),size(D,1),1));
end
