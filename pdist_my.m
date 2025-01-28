function [dist] = pdist_my(X,dist_f)
dist = [];
for i = 1:size(X,1)-1
    dist_ = dist_f(X(i,:,:), X(i+1:end, :,:));
    dist = [dist, dist_]; 
end
end