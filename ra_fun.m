function [ra] = ra_fun(z,tke,z0,HF_option)
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
end
end

