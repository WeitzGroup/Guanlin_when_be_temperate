% model parameters
% Fixed parameters from TS-1 of paper
p.r1_B = 1.4; % growth rate of infected cells
p.K_B = 1e9; % carrying capacity
p.delta_B = 1; % fidelity of vertical transmission
p.beta_B = 200; % burst size
p.m_B = 0.75; % dilution rate
p.b_B = 5e-2; % probability of fusion after adsorption, original paper ~ 1e-2
p.phi_B = 1e-8; % adsorption constant
model_option = 4;

p_set = 0:0.02:1; gamma_set = 1e-2:1e-2:1; % strategy space % Fixed parameters

% R_ver(S):
% direct benefits (R_ver(S) > 1 uniformly) 
% indirect benfits (R_ver(S) < 1 uniformly) 
r1_set_B = [0.2, 3]; 

% Susceptibel host S* is varying in maximum growth rate r
r_set_B = [1+1e-4:1e-4:1+1e-2, 1+2e-2:1e-2:1+1e-1, 1 + 1e-1:1e-1:2, 3:1:3].*p.m_B;
S_set_B = p.K_B.*(1 - p.m_B./r_set_B);

% Rhor and Rver 
Rhor_B = (p.b_B*p.beta_B*p.phi_B).*S_set_B./(p.phi_B.*S_set_B + p.m_B);
Rver_B_ind = (r1_set_B(1).*p.delta_B/p.m_B).*(1 - S_set_B./p.K_B);
Rver_B_di = (r1_set_B(2).*p.delta_B/p.m_B).*(1 - S_set_B./p.K_B);

figure(1);% direct benefits
loglog(S_set_B,Rhor_B,'r','LineWidth',3); hold on;
loglog(S_set_B,Rver_B_di,'b','LineWidth',3); hold on;

% intermediate strategy, p = 0.5, \gamma = 0.1
R0_inter_ind = zeros(length(S_set_B),1); R0_inter_di = zeros(length(S_set_B),1);
pp = 0.5; gg = 0.1; 
for i = 1:length(S_set_B)
    p.r1_B = r1_set_B(1);
    R0_inter_ind(i) = R0_intermediate(p, S_set_B(i), pp, gg, model_option);    
    p.r1_B = r1_set_B(2);
    R0_inter_di(i) = R0_intermediate(p, S_set_B(i), pp, gg, model_option);    
end

figure(1);
loglog(S_set_B,R0_inter_di,'g','LineWidth',3); hold on;
% threhold curve, R0 = 1
loglog(S_set_B,ones(length(S_set_B),1),'k--','LineWidth',2); hold on;


figure(2);% indirect benefits
loglog(S_set_B,Rhor_B,'r','LineWidth',3); hold on;
loglog(S_set_B,Rver_B_ind,'b','LineWidth',3); hold on;
loglog(S_set_B,R0_inter_ind,'g','LineWidth',3); hold on;
% threhold curve, R0 = 1
loglog(S_set_B,ones(length(S_set_B),1),'k--','LineWidth',2); hold on;


figure(1);
axis square; ax = gca; 
xlim([min(S_set_B) max(S_set_B)]); ylim([1e-2 1e2]); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex'); % set Xtick/Ytick in latex interp
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Basic reproduction number, $\mathcal{R}_{0}$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$, direct benefits of vertical transmission',...
    '$\mathcal{R}_{0}$, intermediate strategy'},...
    'Position',[0.52 0.75 0.05 0.05],...
    'Interpreter','latex','FontSize',13); legend boxoff; 
ff = figure(1);
ff.Units = 'inches';
Width = 15; Height = 15;
ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];
Name = 'Berg_model_resource_implicit_lysogeny_direct_benefits';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');



figure(2);
axis square; ax = gca; 
xlim([min(S_set_B) max(S_set_B)]); ylim([1e-2 1e2]); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex'); % set Xtick/Ytick in latex interp
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Basic reproduction number, $\mathcal{R}_{0}$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$, indirect benefits of vertical transmission',...
    '$\mathcal{R}_{0}$, intermediate strategy'},...
    'Position',[0.52 0.75 0.05 0.05],...
    'Interpreter','latex','FontSize',13); legend boxoff; 
% output figures in eps 
ff2 = figure(2);
ff2.Units = 'inches';
Width = 15; Height = 15;
ff2.PaperSize = [Width, Height];
ff2.PaperPosition = [0 0 Width, Height];
ff2.Position = [0 0 Width, Height];
Name2 = 'Berg_model_resource_implicit_lysogeny_indirect_benefits';
print(ff2, Name2, '-depsc2','-r600');
print(ff2, Name2, '-dpdf','-r600');

