function[D] = iter_dict(Y,D,SC)
global option

Ndict = size(D,2);
%D = rand(784,400);
%normdata(D);
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
batch_sz = 256;
    thr = 10;maxIter = 50;momentum=0.5;
    iter = 1;lr = 0.001;converge=0;
     while (~converge)
        iter = iter+1;
        rd = randi([1,size(SC,2)],batch_sz,1);
        SCtmp = SC(:,rd);Ytmp = Y(:,rd);
        grads = 2*(Ytmp-D*SCtmp)*SCtmp';
        if(iter ==2)
            gradsp = grads;
        end
        g = grads*momentum + gradsp*(1-momentum);
        s = D+lr*g;
        D = s./repmat( sqrt(sum(s.*s,1)),size(D,1),1);
        disp(norm(grads));
        grad(iter-1) = norm(grads);
        if( norm(grads)<thr || iter >maxIter||lr < 0.00002)
            converge = 1;
        end
        if (mod(iter,50) == 0)
            if (grad(iter-1)/grad(iter-48)>.8)
                lr = lr/10;
            end
        end
        gradsp = grads;
    end
end
