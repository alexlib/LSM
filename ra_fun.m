function [ra] = ra_fun(z,tke,z0,HF_option,U,zeta)
%compute the SHF areodymic resistance
switch HF_option
    case 'Shao'
        kappa=0.4; %Von Karman constant
        C_k=0.15; %empirical parameter ~0.15
        Pr=0.3; %Prandlt Number from Shao et al. 13
        
        %calculate the Subgrid eddy diffusivity
        K_sg=C_k*(sqrt(tke)/kappa);
        
        %Calculate the subgrid eddy diffusivity for a scalar
        K_hsg=K_sg*Pr^(-1);
        
        %calculate the subgrid areodynamic resistance
        ra =(z/K_hsg)*log(z/z0);
    case 'MOST'
        k=0.4;  %von Karman constant
        gamma1 = mean([28,20.3,19]); %Garratt [1992] Appendix 4
        gamma2 = mean([14,11.6]);       %Garratt [1992] Appendix 4
        beta1M = mean([4,6.0]);         %Garratt [1992] Appendix 4
        beta1H = 7.8;                  %Garratt [1992] Appendix 4
        zT = z0/10;                          %Garratt [1992] Sec 3.3.2 pg 54
        d=0;
        if zeta < 0
            x = (1-gamma1*zeta)^(1/4);
            y = (1-gamma2*zeta)^(1/2);
            psi_M = 2*log((1+x)/2) + log((1+x^2)/2)-2*atan(x)+pi/2;
            psi_H = 2*log((1+y)/2);
            
        else
            psi_M = -beta1M*zeta;
            psi_H = -beta1H*zeta;
            
        end
        
        CH=(k^2)/((log((z-d)/z0)-psi_M)*(log((z-d)/zT)-psi_H));  %CHANGED###aerodynamic transfer coefficient for non-neutral according to MOST in Garrett [1992] "The atmospheric boundary layer"
        ra=1/(CH*U);                  %aerodynamic resistance [s/m]
end
end

