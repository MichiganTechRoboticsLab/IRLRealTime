function [ idx ] = findclosept( pt, arr, mindis )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
idx = -1;

mag = @(x) sqrt(x(1)^2 + x(2)^2);

for ii=1:size(arr,1)
    ret = mag([pt(1) - arr(ii,1), pt(2) - arr(ii,2)]);
    if ret <= mindis
        idx = ii;
        break;
        
    end
end

end

