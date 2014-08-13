function [m, b, R2] = IRLReg(x,y)

        %This is code with finds a regression between the x and y.  It also
        %returns the slope and intercept, plus the R^2 from the regression.

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