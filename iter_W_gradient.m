function [W] = iter_W_gradient(L,P,SC)
global option
Nclassifier = length(P);
batch_sz = 256;maxIter =500;
momentum = 0.5;thr = 0.1;
for i =1:Nclassifier
    M = P{i}*SC;
    converge = 0;
    iter = 1;
    W{i} = rand(size(L,1),size(M,1));
    lr=0.02; 
    while (~converge)
        iter = iter+1;
        rd = randi([1,size(SC,2)],batch_sz,1);
        Mtmp = M(:,rd);Ltmp = L(:,rd);
        grads = (W{i}*Mtmp - Ltmp)*Mtmp'; disp(norm(grads));
        grad(iter-1) = norm(grads);
        if( norm(grads)<thr || iter >maxIter||lr < 0.0002)
            converge = 1;
        end
        if(iter ==2)
            gradsp = grads;
        end
        if (mod(iter,50) == 0)
            if (grad(iter-1)/grad(iter-48)>.8)
                lr = lr/10;
            end
        end
        g = grads*momentum + gradsp*(1-momentum);
        W{i} = W{i} - lr*g;
        gradsp = grads;
    end
        
    
end
end