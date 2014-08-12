function [ b ] = arr2num( a, numBytes )

    b = [];
    
    a = dec2hex(a);
    
    for ii=1:numBytes
        b = [a(ii,:) b];
    end
    
    b=hex2num(b);
end