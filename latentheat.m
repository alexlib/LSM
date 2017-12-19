function [LH] = latentheat(BR,H)
%Latent heat computed from bowen ratio
LH = BR^(-1)*H;
end

