function[trainsplit,validsplit,testsplit] = load_data_scene()
    for i = 0:14
        data = load(['/scene15/data',num2str(i),'mat.mat']);
        data.data = data.data./(repmat(sqrt(sum(data.data.*data.data,1)),size(data.data,1),1));
        trainsplit{i+1} = data.data(:,1:100);
        validsplit{i+1} = data.data(:,101:150);
        testsplit{i+1} = data.data(:,151:end);
    end
end