function[D] = iter_dict(Y,D,SC)
global option
option.u = 100;
Ndict = size(D,2);
% for i = 1:Ndict
%     grads = 100*ones(size(Y,1),1);
%     iter = 0;
%     while iter<4&&norm(grads)>1
% 
%         grads = -2*(Y-D*SC)*SC(i,:)';
%         s = D(:,i)-1/option.u*grads;
%         D(:,i) = s/norm(s,2);
%         iter = iter+1;
%     end
% end
% end

    grads = ones(size(D));
    iter = 0;
    while iter<5
        grads = 2*(Y-D*SC)*SC';
        s = D+1/option.u/20*grads;
        D = s./repmat( sqrt(sum(s.*s,1)),size(D,1),1);
        iter = iter+1;
        disp(norm(grads));
    end
end
