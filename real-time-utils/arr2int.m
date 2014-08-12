function [ b ] = arr2int( a, numBytes )

    b = 0;

    for ii=1:numBytes
        b = b+a(ii)*2^((ii-1)*8);
    end

end