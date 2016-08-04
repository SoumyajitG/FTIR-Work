function [trainsplit,validsplit,testsplit] = load_data_FTIR()
d = 'E:\SCET\FTIR Sample Data';
directory = dir(d);
tmp = 1;
for i = 3 : length(directory)
        name = directory(i).name;
        if name(end-1)=='a' &&name(1) ~='b'&&name(1) ~='t'
            data = importdata(name);
            data = reshape(data,[],1506);
           trainsplit{tmp} = data(20:min(180,size(data,1)),:)';
           validsplit{tmp} = data(1:20,:)';
           %testsplit{tmp} = data(end-50:end,:)';
           testsplit{tmp} = trainsplit{tmp}(:,1:50); 
          tmp = tmp +1;
        end
end


end
