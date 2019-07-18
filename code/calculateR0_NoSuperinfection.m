% p: parameters
% X: environmental data (resource explicit)
% strain : resident or mutant ?
function R_0 = calculateR0_NoSuperinfection(p, X, strategy_pair)

% check if the environment dies out
if X(2) < 1e2
    Sr = 0;
else
    Sr = X(2);
end

if X(1) < 1e-3
    Rr = 0;
else
    Rr = X(1);
end

% (scaled) strategy 
pp = strategy_pair(1);
gg_temp = strategy_pair(2)/(strategy_pair(2) + p.dl);

% define two reproduction numbers (independent of strategy) 
R_hor = p.beta*p.eta*p.phi*Sr*p.alpha/...
             ((p.eta + p.di)*(p.phi*Sr + p.m)*(p.alpha + p.de));
psi = p.b0*Rr/(Rr + p.Rmod);         
R_ver = psi/p.dl;

% construct NGM from reproduction number from loops, P1, P2, P3 
P1 = (1 - pp)*R_hor; P2 = (1 - gg_temp)*R_ver; P3 = R_hor*pp*gg_temp;

% compute fitness
R_0 = 0.5*((P1 + P2 + P3) + ...
     sqrt(P1^2 + P2^2 + P3^2 + 2*P1*P3 + 2*P2*P3 - 2*P1*P2));
                
end