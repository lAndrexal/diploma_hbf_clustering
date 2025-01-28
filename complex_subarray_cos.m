function distance = complex_subarray_cos(XI,XJ, subarray_inds)
XI = XI.' / norm(XI.');
for i = 1:size(XJ,1)
    XJ(i,:) = XJ(i,:) / norm(XJ(i,:));
end
XJ = XJ.';

% distance = 1 - XI'*XJ;
proj = 0;
for i = 1:size(subarray_inds,2)
    inds = subarray_inds(:,i);
    proj = proj + abs(XI(inds)'*XJ(inds,:));
end
% proj = abs(XI'*XJ);

distance = 1 - proj;
end