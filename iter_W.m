function [W] = iter_W(L,P,SC)
global option
Nclassifier = length(P);
for i = 1:Nclassifier
    M = P{i}*SC;
    W{i} = L*M'*inv(M*M'+option.gamma/option.beta*eye(size(M,1)));
end
end