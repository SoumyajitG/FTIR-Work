function [D,W,SC] = joint_train(D,L,train,P)
global option
batch_sz = 100;


Nclassifier = length(P);
maxIter =30;
momentum = 0.5;thr = 0.05;
converge = 0;lrD = 0.01; 
iter = 1;
while (~converge)
    rd =  randi([1,size(train,2)],batch_sz,1);
    traindata = train(:,rd); iter = iter+1;
    trainlabel = L(:,rd);
    %% updata SC
    for i = 1:batch_sz
        SC(:,i) = OMP(D,traindata(:,i),option.T);
    end
    
    %% updata D
    gradsdict = 2*(traindata-D*SC)*SC';
    if(iter ==2)
        gradspdict = gradsdict;
    end
    gdict = gradsdict*momentum + gradspdict*(1-momentum);
    s = D+lrD*gdict;
    D = s./repmat( sqrt(sum(s.*s,1)),size(D,1),1);
    disp(['iter:',num2str(iter),'gradients of dict:',num2str(norm(gradsdict))]);
    graddict(iter-1) = norm(gradsdict);
    if( norm(gradsdict)<thr || iter >maxIter||lrD < 0.00002)
        converge = 1;
    end
    if (mod(iter,50) == 0)
        if (graddict(iter-1)/graddict(iter-48)>.8)
            lrD = lrD/10;
        end
    end
    gradspdict = gradsdict;
    
    if iter ==2
         [W] = iter_W(trainlabel,P,SC);
    
    
    for i = 1:1+option.H
         lrW(i) = 0.01;
    end
    end
    %% updata W
    for i =1:Nclassifier
       M = P{i}*SC;
       gradsW{i} = (W{i}*M- trainlabel)*M';
       disp(['gradients of W',num2str(norm(gradsW{i}))]);
       gradW{i}(iter-1) = norm(gradsW{i});
       if(gradW{i}(iter-1))>1000
           W{i} = trainlabel*M'*inv(M*M'+option.gamma/option.beta*eye(size(M,1)));
           disp('restart');
           lrW(i) = lrW(i)/5;
       end
       if( norm(gradsW{i})<thr || iter >maxIter||lrW(i)< 0.0002)
           converge = 1;
       end
       if(iter ==2)
           gradspW{i} = gradsW{i};
       end
       gW =   gradsW{i}*momentum + gradspW{i}*(1-momentum);
         W{i} = W{i} - lrW(i)*gW;
        gradspW{i} = gradsW{i};
        if (mod(iter,50) == 0)
            if (gradW{i}(iter-1)/gradW{i}(iter-48)>.8)
                lrW(i) = lrW(i)/10;
            end
        end
    end
end
end
