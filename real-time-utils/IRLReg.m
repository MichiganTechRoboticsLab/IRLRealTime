function [m, b, R2] = IRLReg(x,y)

        p = polyfit(x,y,1);
        
        yfit = polyval(p,x,y);
        
        m = p(1);
        b = p(2);
        
        yfit = m*x+b;
        
        yresid=y-yfit;
        
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1)*var(y);
        
        R2 = 1-SSresid/SStotal;

end