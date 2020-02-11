% Input: environmental condition - S* or [R*, S*] 
%        prob_lyso - probability of lysogeny, 0 < p < 1
%        induction rate - \gamma > 0
%        p - fixed model parameters
%        model_option - 4 models, 1,2,3,4

function R0 = R0_intermediate(p, env_cond, prob_lyso, induction_rate, model_option)
pp = prob_lyso; gg = induction_rate;

switch model_option
    case 1 % resource_impliicit present model
        S = env_cond;
        P1 = ((1 - pp)*p.alpha/(p.alpha + p.d1))*(p.beta*p.eta/(p.eta + p.d3))*...
                    (p.phi*S/(p.phi*S + p.m));
        P2 = p.r1*(1 - S/p.K)/(gg + p.d2);
        P3 = (pp*p.alpha/(p.alpha + p.d1))*(p.beta*p.eta/(p.eta + p.d3))*...
                    (p.phi*S/(p.phi*S + p.m))*(gg/(gg + p.d2));
        R0 = 0.5*((P1 + P2 + P3) + ...
                    sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
    case 2 % resource_explicit present model
        R =  env_cond(1); S = env_cond(2); % two environmental factors
        psi_mod = p.b0*R/(R + p.Rmod);
        Rhor = p.phi*S*p.beta*p.eta*p.alpha/...
           ((p.alpha + p.de)*(p.phi*S + p.m)*(p.eta + p.di));
        Rver = psi_mod/p.dl;
        gg_temp = gg/(gg + p.de);
        % reproduction number from loops, P1, P2, P3 
        P1 = (1 - pp)*Rhor; P2 = (1 - gg_temp)*Rver; P3 = Rhor*pp*gg_temp;
        R0 = 0.5*((P1 + P2 + P3) + ...
             sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
        if pp == 0
            R0 = Rhor;
        end
    case 3 % special case: resource_explicit present model, no super immunity lysogen invasion
        R =  env_cond(1); S = env_cond(2); Vr = env_cond(3); % three environmental factors
        psi_mod = p.b0*R/(R + p.Rmod);
        Rhor = p.phi*S*p.beta*p.eta*p.alpha/...
           ((p.alpha + p.de)*(p.phi*S + p.m)*(p.eta + p.di));
        Rver = psi_mod/(p.dl + p.phi*Vr*(1 - p.eps));
        gg_temp = gg/(gg + p.de);
        % reproduction number from loops, P1, P2, P3 
        P1 = (1 - pp)*Rhor; P2 = (1 - gg_temp)*Rver; P3 = Rhor*pp*gg_temp;
        R0 = 0.5*((P1 + P2 + P3) + ...
             sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
        if pp == 0
            R0 = Rhor;
        end
        
    case 4 % Berngruber model
        S = env_cond;
        Rhor_B = (p.b_B*p.beta_B*p.phi_B)*S/(p.phi_B*S + p.m_B);
        Rver_B = (p.r1_B*p.delta_B/p.m_B)*(1 - S/p.K_B);
        gg_temp = gg/(gg + p.m_B);
        % reproduction number from loops, P1, P2, P3 
        P1 = (1 - pp)*Rhor_B; P2 = (1 - gg_temp)*Rver_B; P3 = Rhor_B*pp*gg_temp;
        R0 = 0.5*((P1 + P2 + P3) + ...
                sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
        
        
    case 5 % S&Levin model
       R =  env_cond(1); S = env_cond(2); % two environmental factors
       psi_mod = p.b_SL*R/(R + p.Rmod_SL);
       Rhor_SL = p.phi_SL*S*p.beta_SL/p.rho_SL; Rver_SL = psi_mod/(p.rho_SL + p.tau_SL);
       gg_temp = gg/(gg + p.rho_SL + p.tau_SL);
       % reproduction number from loops, P1, P2, P3 
       P1 = (1 - pp)*Rhor_SL; P2 = (1 - gg_temp)*Rver_SL; P3 = Rhor_SL*pp*gg_temp;
       R0 = 0.5*((P1 + P2 + P3) + ...
            sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
end




end