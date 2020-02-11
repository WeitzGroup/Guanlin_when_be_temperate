% Resource-explicit nonlinear dynamics of temperate phage
% No superinfection, i.e., viruses are only absorbed into susceptible hosts

function dxdt = ode_resource_explicit_NoSuperinfection(t,x,p)

% population states
R = x(1); S = x(2); E = x(3); L = x(4); I = x(5); V = x(6);
% total hosts
N = S + E + L + I;
% resources comsumption and growth rates
psi = p.b0*R/(R + p.Rmod); % Monod growth function
bs = (1 - p.alphaS)*psi; % growth rate of S
bl = psi; % growth rate of L
f_uptake = p.e*(bl*L + bs*S); % resources consumption

% Nonlinear ODE system - strategy : (p.pb, p.gamma)
drdt = p.J - p.rho*R - f_uptake;
dsdt = bs*S - p.phi*S*V  - p.ds*S;
dedt = p.phi*S*V - p.alpha*E - p.de*E;
dldt = bl*L + p.pb*p.alpha*E - p.gamma*L - p.dl*L;
didt = (1 - p.pb)*p.alpha*E - p.eta*I + p.gamma*L - p.di*I;
dvdt = p.beta*p.eta*I - p.phi*S*V - p.m*V;

dxdt = [drdt;dsdt;dedt;dldt;didt;dvdt];
end