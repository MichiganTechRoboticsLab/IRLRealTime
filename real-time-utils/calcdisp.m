function [ t,theta ] = calcdisp( a,b )

%setup coefficient matrix
M=zeros(2*size(a,1),4);
%define coefficient matrix

for kk=1:2:size(a,1)*2
    
    ind = round(kk/2);
    
    M(kk,:)   = [a(ind,1) a(ind,2)  1 0];
    M(kk+1,:) = [a(ind,2) -a(ind,1) 0 1];
    
    N(kk)   = b(ind,1);
    N(kk+1)   = b(ind,2);
    
end

%Computer least squares of data
X=M\N';

%back-calculate angle
theta=[acos(X(1)) asin(X(2))];
theta=real(mean(theta));
t = X(3:4);

end

