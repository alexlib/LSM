function [lambda] = latentheat(T)
%Takes temperature [C] and returns value of latent heat of vaporization [J/g]
  lambda = (2.501 - 0.0024 * T)*1000;
  
end

