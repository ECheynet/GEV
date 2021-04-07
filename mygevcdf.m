function [Fu] = mygevcdf(U,k,sigma,mu)
% [Fu] = mygevcdf(U,k,sigma,mu) estimates the Generalized Extreme Value (GEV) 
% cumulative distribution function (cdf) of a vector U, made of measured 
% extrema.
% 
% Input:
% U: Vector of measured extreme values
% k: shape parameter
% sigma: scale parameter
% mu: location parameter
% 
% Output:
% Fu: GEV cdf, with same size as U
% 
% Example:
% 
% K = [0];
% a = 1;
% u = 1;
% U = linspace(0,2,100);
% F = mygevcdf(U,K,a,u);
% 
% 
% 
% Author: E. Cheynet. University of Stavanger, Norway. Last update:
% 31/12/2016
% 
% See also EVCDF, GEVFIT, GEVINV, GEVLIKE, GEVPDF, GEVRND, GEVSTAT, CDF.
% 

z = (U-mu)/sigma;
if k ~=0,
    Fu = exp(-(1+k.*z).^(-1/k));
else
    Fu = exp(-exp(-z));
end

end

