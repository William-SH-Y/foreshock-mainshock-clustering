function [coef, model] = onestep_linfit(x, y) %coef are the coefficients of the linear fit; model is the line.
     [temp_coef, S] = polyfit(x, y, 1); %coefficients     
     temp_r2 = 1 - (S.normr/norm( y - mean( y )))^2; %r^2
     if isnan(temp_r2) == 1
         temp_r2 = 0; %some r2 are not NaN since the fitting data and the model coefficients are zero
     end
     coef = [temp_coef, temp_r2]; %The third number is r2
     
     model = coef(1)*x + coef(2);
    
end