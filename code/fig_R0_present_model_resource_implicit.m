% model parameters
p.r1 = 1.2;  % max growth rate of lysogenic cells
p.d = 0.2; % death rate of susceptible cells
p.d1 = 0.2; % death rate of infected cells
p.d2 = 0.2; % death rate of lysogenic cells
p.d3 = 0.2; % death rate of lytic cells
p.m = 1/24; % decay rate of phage
p.K = 2e8; % carrying capacity 
p.phi = 3.4e-10;  % adsorption rate
p.alpha = 2; % transit decision rate
p.eta = 1; % lysis rate
p.beta = 50; % burst rate

% viral strategy space
p_set = 0:0.02:1; gamma_set = 1e-2:1e-2:1; 

% range of susceptible host (vary with free parameter, growth rate of S)
r_set = [1+5e-4:2e-4:1+1e-2, 1+2e-2:1e-2:1+1e-1, 1+2e-1:1e-1:5].*p.d;
S_set = p.K.*(1 - p.d./r_set);

% R_hor(S) 
Rhor = (p.beta*p.eta*p.phi*p.alpha).*S_set./...
    ((p.eta + p.d3).*(p.phi.*S_set + p.m).*(p.alpha + p.d1));

% R_ver(S):
% direct benefits (R_ver(S) > 1 uniformly) 
% indirect benfits (R_ver(S) < 1 uniformly) 
r1_set = [0.1, 1.1]; 
Rver_ind = (r1_set(1)/p.d2).*(1 - S_set./p.K);
Rver_di = (r1_set(2)/p.d2).*(1 - S_set./p.K);
figure(1);% direct benefits
loglog(S_set,Rhor,'r','LineWidth',3); hold on;
loglog(S_set,Rver_di,'b','LineWidth',3); hold on;

%{
tmpt1 = text(10^5.1,10^-0.55,{'Indirect benefits';'of lysogeny'});
set(tmpt1,'fontsize',12,'interpreter','latex');
tmpt1 = text(10^5.1,10^0.6,{'Direct benefits of lysogeny'});
set(tmpt1,'fontsize',12,'interpreter','latex');
%}

% intermediate strategy, p = 0.5, \gamma = 0.1
R0_inter_ind = zeros(length(S_set),1); R0_inter_di = zeros(length(S_set),1);
pp = 0.5; gg = 0.1; model_opt = 1;
for i = 1:length(S_set)
    p.r1 = r1_set(1);
    R0_inter_ind(i) = R0_intermediate(p, S_set(i), pp, gg, model_opt);    
    p.r1 = r1_set(2);
    R0_inter_di(i) = R0_intermediate(p, S_set(i), pp, gg, model_opt);    
end
figure(1);
loglog(S_set,R0_inter_di,'g','LineWidth',3); hold on;
% threhold curve, R0 = 1
loglog(S_set,ones(length(S_set),1),'k--','LineWidth',2); hold on;


figure(2);% indirect benefits
loglog(S_set,Rhor,'r','LineWidth',3); hold on;
loglog(S_set,Rver_ind,'b','LineWidth',3); hold on;
loglog(S_set,R0_inter_ind,'g','LineWidth',3); hold on;
% threhold curve, R0 = 1
loglog(S_set,ones(length(S_set),1),'k--','LineWidth',2); hold on;

figure(1);
axis square; ax = gca; 
xlim([min(S_set) max(S_set)]); ylim([1e-2 1e2]); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex'); % set Xtick/Ytick in latex interp
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Basic reproduction number, $\mathcal{R}_{0}$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$, direct benefits of vertical transmission',...
    '$\mathcal{R}_{0}$, intermediate strategy'},...
    'Position',[0.50 0.75 0.05 0.05],...
    'Interpreter','latex','FontSize',13); legend boxoff; 

% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 15; Height = 15;
ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];
Name = 'present_model_resource_implicit_lysogeny_direct_benefits';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');

figure(2);
axis square; ax = gca; 
xlim([min(S_set) max(S_set)]); ylim([1e-2 1e2]); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex'); % set Xtick/Ytick in latex interp
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Basic reproduction number, $\mathcal{R}_{0}$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$, indirect benefits of vertical transmission',...
    '$\mathcal{R}_{0}$, intermediate strategy'},...
    'Position',[0.50 0.75 0.05 0.05],...
    'Interpreter','latex','FontSize',13); legend boxoff; 
% output figures in eps 
ff2 = figure(2);
ff2.Units = 'inches';
Width = 15; Height = 15;
ff2.PaperSize = [Width, Height];
ff2.PaperPosition = [0 0 Width, Height];
ff2.Position = [0 0 Width, Height];
Name2 = 'present_model_resource_implicit_lysogeny_indirect_benefits';
print(ff2, Name2, '-depsc2','-r600');
print(ff2, Name2, '-dpdf','-r600');
