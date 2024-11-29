function [paramH, yextra] = hill_eq(paramH_0,x,y,xextra)

%xextra:optional
if nargin==3
    xextra=x;
end

fun= @(param,x) (param(1)) ./ (1+(param(2)./x).^param(3)); %hill equation

paramH = lsqcurvefit(fun,paramH_0,x,y);

yextra = fun(paramH,xextra);
end