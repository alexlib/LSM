function [saturated] = satvap(T,H,BR)
%Calculate the saturated vapor pressure
K = T+273.15; %[K]
lambda = 1000*latentheat(T);
saturated = 611*exp((lambda/461)*((1/273.15)-(1/K)));
end

