function [ idx, closest ] = findclosepts( arr, pt, k, mindis )
idx = [];

mag = @(x) sqrt(x(1)^2 + x(2)^2);

for ii=1:size(arr,1)
    
    if ii==k
        continue
    end
    
    ret = mag([pt(1) - arr(ii,1), pt(2) - arr(ii,2)]);
    if ret <= mindis
        idx = [idx ii];
    end
end

ptdiffs = [];
for ii=1:length(idx)
    
    sub = arr(idx(ii),:) - pt;
    ptdiffs = [ptdiffs ; mag(sub)];
    
end
[~,I] = min(ptdiffs);
closest = idx(I);

end

