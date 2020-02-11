% Fixed parameters in Stewart & Levins (1984)
p.beta_SL = 100; % burst size
p.phi_SL = 1e-9; % absorption rate
p.tau_SL = 1e-3; % loss of prophage
p.Rmod_SL = 4; % half saturation
p.rho_SL = 0.2; % dilution rate 
p.e_SL = 5e-7; % conversion coefficient
p.b_SL = 0.7; % maximum growth rate
p.gamma_min = 1e-2;

% S* and R* are varying in selection coeffs alpha_s and C
alphaS_set = 0.1:1e-2:(1 - p.rho_SL/p.b_SL) - 0.1; 
% theta_set = 1e-2:1e-2:10; % auxillary free var, theta = C - R*(alpha_s)
theta_set = 2e-2:1e-2:50;
data_C = zeros(length(alphaS_set),length(theta_set));
% range of environments
R_set = p.Rmod_SL.*p.rho_SL./((1 - alphaS_set).*p.b_SL - p.rho_SL);
S_set = theta_set./p.e_SL;

tau_SL_set = [1e-3 1];

% intermediate lysogeny probabilty and induction rate
pp = 0.5; gg = 1e-1; model_option = 5;
% Save horizontal and vertical transmission fitness at virus-free
data_Rver_ind = zeros(length(alphaS_set),length(theta_set));
data_Rver_di = zeros(length(alphaS_set),length(theta_set));
data_Rhor = zeros(length(alphaS_set),length(theta_set));
% Save intermediate transmission fitness at virus-free
data_R0_inter_ind = zeros(length(alphaS_set),length(theta_set));
data_R0_inter_di = zeros(length(alphaS_set),length(theta_set));

for i = 1:length(alphaS_set)
    for j = 1:length(theta_set)
        S_SL = S_set(j); R_SL = R_set(i);
        % assign two model parameters
        psi_mod = p.b_SL*R_SL/(R_SL + p.Rmod_SL);
        % theta = C - R*(alpha_s)
        data_C(i,j) = R_SL*alphaS_set(i) + theta_set(j);
        % two exclusive strategy fitness
        data_Rhor(i,j) = p.phi_SL*S_SL*p.beta_SL/p.rho_SL;
        data_Rver_di(i,j) = psi_mod/(p.rho_SL + tau_SL_set(1));
        data_Rver_ind(i,j) = psi_mod/(p.rho_SL + tau_SL_set(2));
               
        % compute intermediate R0
        p.tau_SL = tau_SL_set(1); % direct benefits
        R0_inter_di = ...
            R0_intermediate(p, [R_SL, S_SL], pp, gg, model_option);
        p.tau_SL = tau_SL_set(2); % indirect benefits
        R0_inter_ind = ...
            R0_intermediate(p, [R_SL, S_SL], pp, gg, model_option);
        
        % save two fitness values
        data_R0_inter_ind(i,j) = R0_inter_ind;
        data_R0_inter_di(i,j) = R0_inter_di;        
    end
end

% define unified color bar for all plots
lower_lim_cbar = min([round(min(min(log10(data_R0_inter_ind))),1), round(min(min(log10(data_R0_inter_di))),1), ...
    round(min(min(log10(data_Rhor))),1), round(min(min(log10(data_Rver_di))),1)]);
upper_lim_cbar = max([round(max(max(log10(data_R0_inter_ind))),1), round(max(max(log10(data_R0_inter_di))),1),...
    round(max(max(log10(data_Rhor))),1), round(max(max(log10(data_Rver_di))),1)]);


%% --------------- Heatmap of horizontal fitness ----------------------------
figure(1);
%heatmap of horizontal fitness
pcolor(S_set,R_set,log10(data_Rhor)); 
h = pcolor(S_set,R_set,log10(data_Rhor)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
shading interp;
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c1 = caxis;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);

% contour for {R0 = 1}
levels =  [0 0];
[~,c] = contour(S_set,R_set,log10(data_Rhor),levels, '--k');
c.LineWidth = 3;
tmpt1 = text(10^7, 5.5, {'$\mathcal{R}^{hor}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt1 = text(1e5, 5.5, {'$\mathcal{R}^{hor} < 1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):1:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Horizontal invasion fitness, $\mathcal{R}^{hor}$',...
    'fontsize',20,'interpreter','latex');



%% --------------- Heatmap of vertical fitness (direct benefits) ----------------------------
figure(2);
pcolor(S_set,R_set,log10(data_Rver_di)); 
h = pcolor(S_set,R_set,log10(data_Rver_di)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
shading interp;
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c2 = caxis;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);
tmpt1 = text(10^6, 6.9, {'$\mathcal{R}^{ver}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):1:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Vertical invasion fitness, $\mathcal{R}^{ver}$',...
    'fontsize',20,'interpreter','latex');




%% --------------- Heatmap of vertical fitness (indirect benefits) ----------------------------
figure(3);
pcolor(S_set,R_set,log10(data_Rver_ind)); 
h = pcolor(S_set,R_set,log10(data_Rver_ind)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c3 = caxis;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);
% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
tmpt1 = text(10^6, 6.9, {'$\mathcal{R}^{ver}<1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):1:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Vertical invasion fitness, $\mathcal{R}^{ver}$',...
    'fontsize',20,'interpreter','latex')


%% --------------- Heatmap of intermediate strategy fitness (w indirect benefits) ----------------------------
figure(4);
pcolor(S_set,R_set,log10(data_R0_inter_ind));
h = pcolor(S_set,R_set,log10(data_R0_inter_ind)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
shading interp;
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c4 = caxis;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);

% contour for {R0 = 1}
levels =  [0 0];
[~,c] = contour(S_set,R_set, log10(data_R0_inter_ind), levels, '--k');
c.LineWidth = 3;
tmpt1 = text(10^7, 5.5, {'$\mathcal{R}_{0}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt1 = text(1e5, 5.5, {'$\mathcal{R}_{0} < 1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');

% lables, axises
set(gca,'XScale','log'); 
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):1:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Invasion fitness with $0 < p < 1$, $\mathcal{R}_{0}$',...
    'fontsize',20,'interpreter','latex');
axis square;




%% --------------- Heatmap of intermediate strategy fitness (w direct benefits) ----------------------------
figure(5);
pcolor(S_set,R_set,log10(data_R0_inter_di)); 
h = pcolor(S_set,R_set,log10(data_R0_inter_di)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
shading interp;
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c5 = caxis;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);

% contour for {R0 = 1}
levels = 1;
[~,c] = contour(S_set,R_set,log10(data_R0_inter_di), [0 0], '--k');
c.LineWidth = 3;
tmpt1 = text(1e7, 9, {'$\mathcal{R}_{0}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt2 = text(1e5, 3.9, {'$\mathcal{R}_{0}<1$'});
set(tmpt2,'fontsize',25,'interpreter','latex');
x = [0.30 0.26]; % start
y = [0.29 0.2]; % end
annotation('arrow',x,y, 'LineWidth',2, 'HeadLength',15, 'HeadWidth',15)

% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):1:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Invasion fitness with $0 < p < 1$, $\mathcal{R}_{0}$',...
    'fontsize',20,'interpreter','latex');


% unifying all the cbars
c_f = [min([c1 c2 c3 c4 c5]), max([c1 c2 c3 c4 c5])]; 
for i = 1:5
    figure(i);
    caxis(c_f);
end


%% --------------- Output files ----------------------------
ff1 = figure(1);
ff1.Units = 'inches';
Width = 17; Height = 17;
ff1.PaperSize = [Width, Height];
ff1.PaperPosition = [0 0 Width, Height];
ff1.Position = [0 0 Width, Height];
Name1 = 'SB_model_resource_explicit_horizontal_fitness';
print(ff1, Name1, '-depsc2','-r600');
print(ff1, Name1, '-dpdf','-r600');


% Output files
ff2 = figure(2);
ff2.Units = 'inches';
Width = 17; Height = 17;
ff2.PaperSize = [Width, Height];
ff2.PaperPosition = [0 0 Width, Height];
ff2.Position = [0 0 Width, Height];
Name2 = 'SB_model_resource_explicit_vertical_fitness_direct_benefits';
print(ff2, Name2, '-depsc2','-r600');
print(ff2, Name2, '-dpdf','-r600');

% Output files
ff3 = figure(3);
ff3.Units = 'inches';
Width = 17; Height = 17;
ff3.PaperSize = [Width, Height];
ff3.PaperPosition = [0 0 Width, Height];
ff3.Position = [0 0 Width, Height];
Name3 = 'SB_model_resource_explicit_lysogeny_indirect_benefits';
print(ff3, Name3, '-depsc2','-r600');
print(ff3, Name3, '-dpdf','-r600');

% Output files
ff4 = figure(4);
ff4.Units = 'inches';
Width = 17; Height = 17;
ff4.PaperSize = [Width, Height];
ff4.PaperPosition = [0 0 Width, Height];
ff4.Position = [0 0 Width, Height];
Name4 = 'SB_model_resource_explicit_intermed_fitness_indirect_benefits';
print(ff4, Name4, '-depsc2','-r600');
print(ff4, Name4, '-dpdf','-r600');

ff5 = figure(5);
ff5.Units = 'inches';
Width = 17; Height = 17;
ff5.PaperSize = [Width, Height];
ff5.PaperPosition = [0 0 Width, Height];
ff5.Position = [0 0 Width, Height];
Name5 = 'SB_model_resource_explicit_intermed_fitness_direct_benefits';
print(ff5, Name5, '-depsc2','-r600');
print(ff5, Name5, '-dpdf','-r600');

