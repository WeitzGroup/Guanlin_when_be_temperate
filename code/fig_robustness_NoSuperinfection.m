% Model parameters
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
p.m = 0.4;

% read data in, evolution robustness data
evo_robustness_data = importdata('data_robustness_evo_robustness_resource_explicit_NoSuperinfection.txt');
theta_set = [0.01:0.01:0.1, 0.2:0.05:1, 2:0.5:15]; % augumented parameter range
alphaS_set = -1:0.02:0.4; % alphaS < 0, S grows faster than L; alphaS > 0, ...

% Range of Susceptible hosts and resources
S_set = theta_set./(p.e*p.ds);
R_set = p.Rmod*p.ds./((1 - alphaS_set).*p.b0 - p.ds);
R_hor_set = p.beta*p.eta*p.phi.*S_set.*p.alpha./...
             ((p.eta + p.di)*(p.phi*S_set + p.m)*(p.alpha + p.de));
R_c = p.Rmod./(p.b0./(R_hor_set.*p.dl) - 1);

maph = [ 97/255    151/255    244/255 % light blue
         0.8 0.8 0.8 % gray
         1 0 0 % red
         0 0 1]; % blue   
figure;
pcolor(S_set,R_set,evo_robustness_data); 
h = pcolor(S_set,R_set,evo_robustness_data); hold on;
set(h, 'EdgeColor', 'none'); 
colormap(maph); alpha(h,0.7);

plot(S_set, R_c, 'k', 'LineWidth',4); hold on;
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20);
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(0.5:0.2:max(R_set)); 
xlabel({'Susceptible population, $S^{*} (ml^{-1})$ '}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ylabel({'Resources, $R^{*} (\mu g/ml)$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
     