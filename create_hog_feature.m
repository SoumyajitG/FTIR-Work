function [] = create_hog_feature()
[trainsplit,validsplit,testsplit] = load_data_mnist();
for i = 1:10
    data = [];
    for j = 1:size(trainsplit{i},2)
        im = reshape(trainsplit{i}(:,j),28,28);
        hog = hog_feature_vector (im);
        data = [data,hog'];
        if mod(j,100)==0
            disp(j);
        end
    end
    for j = 1:size(validsplit{i},2)
        im = reshape(validsplit{i}(:,j),28,28);
        hog = hog_feature_vector (im);
        data = [data,hog'];
         if mod(j,100)==0
            disp(j);
        end
    end
    train = data'*255;
    
    save(['train',num2str(i-1),'.mat'],'train');
    data = [];
    for j = 1:size(testsplit{i},2)
        im = reshape(testsplit{i}(:,j),28,28);
        hog = hog_feature_vector (im);
        data = [data,hog'];
         if mod(j,100)==0
            disp(j);
        end
    end
    test = data'*255;
    save(['test',num2str(i-1),'.mat'],'test');
end
end

