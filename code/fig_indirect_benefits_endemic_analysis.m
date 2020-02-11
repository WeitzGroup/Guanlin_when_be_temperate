%% Model parameters
% Note: same as the evolutionary invasion analysis section
% resource parameters
p.Rmod = 5; % half saturation
p.rho = 0.5; % dilution rate 
p.e = 5e-7; % conversion coefficient
p.b0 = 1.2; % maximum growth rate
% explicit infection parameters 
p.eta = 1; % lysis rate
p.beta = 50; % burst rate
p.phi = 3.4e-10;  % adsorption rate
p.alpha = 2; % transit decision rate
% cell traits
p.ds = 0.2; % death rate of susceptible cells
p.de = 0.2; % death rate of infected cells
p.dl = 0.2; % death rate of lysogenic cells
p.di = 0.2; % death rate of lytic cells
p.m = 1/24; % phage decay rate

% (R*, S*) is varied by the parameter space (alphaS, J);
theta_set = [0.01:0.01:0.1, 0.2:0.05:1, 2:0.5:10]; % augumented parameter range
alphaS_set = -1:0.02:0.5; % alphaS < 0, S grows faster than L; alphaS > 0, ...
% Range of Susceptible hosts and resources
S_set = theta_set./(p.e*p.ds);
R_set = p.Rmod*p.ds./((1 - alphaS_set).*p.b0 - p.ds);


dl_set = [0.8, 0.08]; % indirect benefits and direct benefits

equilibrium_id_endemic = importdata('data_equilibrium_id_endemic.txt');
vertical_fitness_endemic = importdata('data_vertical_fitness_endemic.txt');

% classify data 
[row0,col0] = find(vertical_fitness_endemic > 0 & vertical_fitness_endemic < 1);
[row1,col1] = find(vertical_fitness_endemic > 1);
[row2,col2] = find(vertical_fitness_endemic < 0);
inv_id_endemic = vertical_fitness_endemic;

% indirect benefits lysogen cannot invade
for i = 1:length(row0)
    inv_id_endemic(row0(i),col0(i)) = -3;
end
% indirect benefits lysogen invade
for i = 1:length(row1)
    inv_id_endemic(row1(i),col1(i)) = -1;
end
for i = 1:length(row2) % oscillatory
    inv_id_endemic(row2(i),col2(i)) = -2;
end

figure(1); % invasion of lytic viruses endemic state via indirect benefits lysogens
h = pcolor(S_set,R_set,inv_id_endemic); 
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.5); colormap('flag'); 

set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 

tmpt1 = text(10^5.3, 1.5, {'$\mathcal{R}^{hor} < 1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
%tmpt2 = text(10^7.7, 2.2, {'fluctuating';'dynamics'});
%set(tmpt2,'fontsize',15,'interpreter','latex');

xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title({'Invasion of resident lytic virus via a'; 'strain with indirect vertical fitness benefits'},...
    'fontsize',20,'interpreter','latex');



