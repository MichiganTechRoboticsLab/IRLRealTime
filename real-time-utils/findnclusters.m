function [ arr ] = findnclusters( arr, ii, n, dx )
%This function marks an array of points with group numbers.  It is a
%recurisve algorithm, which is basically a depth-first search.

marked = @(x,ii)   x(ii,3)~=0;

if size(arr,2) < 3
    arr=[arr zeros(size(arr,1),1)];
end

if marked(arr,ii)
    return
end

arr(ii,3) = n;
idx = findclosepts(arr,arr(ii,:), ii, dx);
for jj=1:length(idx)
    arr=findnclusters(arr,idx(jj),n,dx);
end

ii = find(arr(:,3) == 0, 1);

if(~isempty(ii))
    arr = findnclusters(arr,ii,n+1, dx);
end

end
