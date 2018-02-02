function [ra] = ra_fun(z,tke,z0,HF_option,u_star,zeta,Pr)
%compute the SHF areodymic resistance
%in the case of bulk, pass u, the mean wind speed instead of u_star
switch HF_option
    case 'Shao'
        kappa=0.4; %Von Karman constant
        C_k=0.15; %empirical parameter ~0.15
        %Prandlt Number from Shao et al. 13
        
        %calculate the Subgrid eddy diffusivity
        K_sg=C_k.*sqrt(tke).*z;
        
        %Calculate the subgrid eddy diffusivity for a scalar
        K_hsg=K_sg.*Pr.^(-1);
        
        %calculate the subgrid areodynamic resistance
        ra =(z./K_hsg).*log(z./z0);
    case 'Bulk'
        k=0.4;  %von Karman constant
        gamma1 = 19; %Garratt [1992] Appendix 4, Hogstrom 1988
        gamma2 = 11.6;       %Garratt [1992] Appendix 4
        beta1M =6;         %Garratt [1992] Appendix 4
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
        ra=1/(CH*u_star);                  %aerodynamic resistance [s/m]
    case 'MOST'
        k=0.4;  %von Karman constant
        gamma1 = 19; %Garratt [1992] Appendix 4, Hogstrom 1988
        gamma2 = 11.6;       %Garratt [1992] Appendix 4
        beta1M =6;         %Garratt [1992] Appendix 4
        beta1H = 7.8;                  %Garratt [1992] Appendix 4
        zT = z0/10;                          %Garratt [1992] Sec 3.3.2 pg 54
        d=0;
        if zeta < 0
            x = (1-gamma1.*zeta).^(1/4);
            y = (1-gamma2.*zeta).^(1/2);
            psi_M = 2.*log((1+x)/2) + log((1+x.^2)./2)-2.*atan(x)+pi./2;
            psi_H = 2.*log((1+y)./2);
            
        else
            psi_M = -beta1M.*zeta;
            psi_H = -beta1H.*zeta;
            
        end
        ra = (1./(k.*u_star)).*(log(z./zT) - psi_H);
        %aerodynamic resistance [s/m]
end
end

