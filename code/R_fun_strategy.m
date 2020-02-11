% Input: parameters (p), environment condition (e.g., low S, high S)
%        strategy space grid, p and gamma
%        model, {1,2,3}, 1 - present model, 2 - B model, 3 - S&L model
% Output: R_matrix is the function of p and gamma given environmental
%         condition.

function R_matrix = R_fun_strategy(p, envir, p_set, gamma_set, model)
R_matrix = zeros(length(p_set),length(gamma_set)); % initialization

switch model
    case 1 % present model
        S = envir;
        for i = 1:length(p_set)
            for j = 1:length(gamma_set)
                pp = p_set(i); gg = gamma_set(j);
                % reproduction number from loops, P1, P2, P3 
                P1 = ((1 - pp)*p.alpha/(p.alpha + p.d1))*(p.beta*p.eta/(p.eta + p.d3))*...
                    (p.phi*S/(p.phi*S + p.m));
                P2 = p.r1*(1 - S/p.K)/(gg + p.d2);
                P3 = (pp*p.alpha/(p.alpha + p.d1))*(p.beta*p.eta/(p.eta + p.d3))*...
                    (p.phi*S/(p.phi*S + p.m))*(gg/(gg + p.d2));
                R0 = 0.5*((P1 + P2 + P3) + ...
                    sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
                R_matrix(i,j) = R0;
            end
        end
        
    case 2 % Berngruber model
        S = envir;
        Rhor_B = (p.b_B*p.beta_B*p.phi_B)*S/(p.phi_B*S + p.m_B);
        Rver_B = (p.r1_B*p.delta_B/p.m_B)*(1 - S/p.K_B);
        for i = 1:length(p_set)
            for j = 1:length(gamma_set)
                pp = p_set(i); gg = gamma_set(j); gg_temp = gg/(gg + p.m_B);
                % reproduction number from loops, P1, P2, P3 
                P1 = (1 - pp)*Rhor_B; P2 = (1 - gg_temp)*Rver_B; P3 = Rhor_B*pp*gg_temp;
                R0 = 0.5*((P1 + P2 + P3) + ...
                    sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
                R_matrix(i,j) = R0;
            end
        end
        
    case 3 % S&Levin model
       R =  envir(1); S = envir(2); % two environmental factors
       psi_mod = p.b_SL*R/(R + p.Rmod_SL);
       Rhor_SL = p.phi_SL*S*p.beta_SL/p.rho_SL; Rver_SL = psi_mod/(p.rho_SL + p.tau_SL);
       for i = 1:length(p_set)
            for j = 1:length(gamma_set)
                pp = p_set(i); gg = gamma_set(j); gg_temp = gg/(gg + p.rho_SL + p.tau_SL);
                % reproduction number from loops, P1, P2, P3 
                P1 = (1 - pp)*Rhor_SL; P2 = (1 - gg_temp)*Rver_SL; P3 = Rhor_SL*pp*gg_temp;
                R0 = 0.5*((P1 + P2 + P3) + ...
                    sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
                R_matrix(i,j) = R0;
            end
       end
       
    case 4 % Present model, resource explicit version
       R =  envir(1); S = envir(2); % two environmental factors
       psi_mod = p.b0*R/(R + p.Rmod);
       Rhor = p.phi*S*p.beta*p.eta*p.alpha/...
           ((p.alpha + p.de)*(p.phi*S + p.m)*(p.eta + p.dl));
       Rver = psi_mod/p.dl;
       for i = 1:length(p_set)
            for j = 1:length(gamma_set)
                pp = p_set(i); gg = gamma_set(j); gg_temp = gg/(gg + p.de);
                % reproduction number from loops, P1, P2, P3 
                P1 = (1 - pp)*Rhor; P2 = (1 - gg_temp)*Rver; P3 = Rhor*pp*gg_temp;
                R0 = 0.5*((P1 + P2 + P3) + ...
                    sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
                R_matrix(i,j) = R0;
            end
       end
       
       
end
end
