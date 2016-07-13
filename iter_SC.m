function[SC]=iter_SC(dict,L,Y,W,P)
% dict  featuredim * Ndict
% W  Ncls * Nfeature
% P Nfeature*Ndict
global option
Nclassifier = length(W);
for Nsample = 1:size(Y,2)
    newY = [Y(:,Nsample);sqrt(option.beta)*repmat(L(:,Nsample),Nclassifier,1)];
    newA = dict;
    for i = 1:Nclassifier
        newA = [newA;sqrt(option.beta)*W{i}*P{i}];
    end
    SC(:,Nsample) = OMP(newA,newY,option.T);
end
end
