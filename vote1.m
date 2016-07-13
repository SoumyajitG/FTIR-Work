function [cls] = vote1(score)
score = sum(score,2);
id = find(score == max(score));
cls = zeros(size(score,1),1);
cls(id(1)) = 1;
end