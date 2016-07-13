function[trainsplit,validsplit,testsplit] = load_data_mnist()
    for i = 0:9
        data = load(['mnist/train',num2str(i),'.mat']);
        trainsplit{i+1} = data.D(1:4000,:)'/255;
        validsplit{i+1} = data.D(4001:end,:)'/255;
        data = load(['mnist/test',num2str(i),'.mat']);
        testsplit{i+1} = data.D'/255;
    end
end