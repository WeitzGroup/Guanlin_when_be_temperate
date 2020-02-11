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

% assume gamma from [1e-2, 1];
p.gamma_min = 1e-2; p.gamma_max = 1;

% (R*, S*) is varied by the parameter space (alphaS, J);
theta_set = [0.01:0.01:0.1, 0.2:0.05:1, 2:0.5:10]; % augumented parameter range
alphaS_set = -1:0.02:0.5; % alphaS < 0, S grows faster than L; alphaS > 0, ...
data_alphaS = zeros(length(alphaS_set),length(theta_set));
data_J = zeros(length(alphaS_set),length(theta_set));
% Range of Susceptible hosts and resources
S_set = theta_set./(p.e*p.ds);
R_set = p.Rmod*p.ds./((1 - alphaS_set).*p.b0 - p.ds);

% R_hor_set(S*), using S* to compute the range of horizonral fitness
R_hor_set = p.beta*p.eta*p.phi.*S_set.*p.alpha./...
             ((p.eta + p.di)*(p.phi*S_set + p.m)*(p.alpha + p.de));
         
% R_ver(R*):
% direct benefits (R_ver(R*) > 1 uniformly), indirect benfits (R_ver(R*) < 1 uniformly) 
% we use decay rate of lysogen, dl, to distinguish two cases
dl_set = [0.5, 0.08]; % indirect benefits and direct benefits
%psi_star_set = p.b0.*R_set./(R_set + p.Rmod);         
%R_ver_ind_set = psi_star./dl_set(1);
%R_ver_di_set = psi_star./dl_set(2);
        
% Save horizontal and vertical transmission fitness at virus-free
data_Rver_ind = zeros(length(alphaS_set),length(theta_set));
data_Rver_di = zeros(length(alphaS_set),length(theta_set));
data_Rhor = zeros(length(alphaS_set),length(theta_set));
% Save intermediate transmission fitness at virus-free
data_R0_inter_ind = zeros(length(alphaS_set),length(theta_set));
data_R0_inter_di = zeros(length(alphaS_set),length(theta_set));

% intermediate lysogeny probabilty and induction rate
pp = 0.5; gg = 1e-1; model_option = 2;

for i = 1:length(alphaS_set)
    for j = 1:length(theta_set)
        % assign two model parameters
        p.alphaS = alphaS_set(i);
        p.J = theta_set(j) + (p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds))*p.rho;
        % Save the parameter pair
        data_alphaS(i,j) = p.alphaS; data_J(i,j) = p.J;
        
        % compute virus-free environmental states
        R_star = p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        S_star = theta_set(j)/(p.e*p.ds);

        % compute two transmission fitness values
        R_hor = p.beta*p.eta*p.phi*S_star*p.alpha/...
             ((p.eta + p.di)*(p.phi*S_star + p.m)*(p.alpha + p.de));
        psi_star = p.b0*R_star/(R_star + p.Rmod);      
        R_ver_ind = psi_star/dl_set(1); R_ver_di = psi_star/dl_set(2);
        
        % compute intermediate R0
        p.dl = dl_set(1); % indirect benefits
        R0_inter_ind = ...
            R0_intermediate(p, [R_star, S_star], pp, gg, model_option);
        p.dl = dl_set(2); % direct benefits
        R0_inter_di = ...
            R0_intermediate(p, [R_star, S_star], pp, gg, model_option);
        
        % save two fitness values
        data_Rver_ind(i,j) = R_ver_ind; data_Rver_di(i,j) = R_ver_di; 
        data_Rhor(i,j) = R_hor;       
        data_R0_inter_ind(i,j) = R0_inter_ind;
        data_R0_inter_di(i,j) = R0_inter_di;        
    end
end

% define unified color bar for all plots
lower_lim_cbar = min([round(min(min(log10(data_R0_inter_ind))),1), round(min(min(log10(data_R0_inter_di))),1), ...
    round(min(min(log10(data_Rhor))),1), round(min(min(log10(data_Rver_di))),1)]);
upper_lim_cbar = max([round(max(max(log10(data_R0_inter_ind))),1), round(max(max(log10(data_R0_inter_di))),1),...
    round(max(max(log10(data_Rhor))),1), round(max(max(log10(data_Rver_di))),1)]);

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
tmpt1 = text(10^7, 1.5, {'$\mathcal{R}^{hor}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt1 = text(2*10^5, 1.5, {'$\mathcal{R}^{hor} < 1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Horizontal invasion fitness, $\mathcal{R}^{hor}$',...
    'fontsize',20,'interpreter','latex');
%{
tmptA = text(7.5e7, 1.28, {'$A$'});
set(tmptA,'fontsize',20,'interpreter','latex');
tmptB = text(5e7, 1.10, {'$B$'});
set(tmptB,'fontsize',20,'interpreter','latex');
tmptC = text(5e7, 0.66, {'$C$'});
set(tmptC,'fontsize',20,'interpreter','latex');
%}

figure(2);
%heatmap of horizontal fitness
pcolor(S_set,R_set,log10(data_Rver_di)); 
h = pcolor(S_set,R_set,log10(data_Rver_di)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
c2 = caxis;
shading interp;
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);
tmpt1 = text(10^6, 1.5, {'$\mathcal{R}^{ver}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');

% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Vertical invasion fitness, $\mathcal{R}^{ver}$',...
    'fontsize',20,'interpreter','latex');

figure(3);
%heatmap of horizontal fitness
pcolor(S_set,R_set,log10(data_Rver_ind)); 
h = pcolor(S_set,R_set,log10(data_Rver_ind)); hold on;
set(h, 'EdgeColor', 'none'); hold on; alpha(h,0.8); colormap('hot'); 
cbar = colorbar; cbary = get(cbar,'ylabel'); cbar.LineWidth = 2;
c3 = caxis;
shading interp;
set(cbar,'TickLabelInterpreter', 'latex');
set(cbary,'string','Log$_{10}$ Fitness','interpreter','latex',...
    'FontName', 'Times New Roman','FontSize',20');
set(cbar, 'Ticks', lower_lim_cbar:0.3:upper_lim_cbar);
tmpt1 = text(10^6, 1.5, {'$\mathcal{R}^{ver}<1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');

% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Vertical invasion fitness, $\mathcal{R}^{ver}$',...
    'fontsize',20,'interpreter','latex')


figure(4);
%heatmap of horizontal fitness
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
tmpt1 = text(10^7, 1.5, {'$\mathcal{R}_{0}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt2 = text(2e5, 1.5, {'$\mathcal{R}_{0}<1$'});
set(tmpt2,'fontsize',25,'interpreter','latex');

% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 
xlabel({'Susceptible population, $S^{*}$ (ml$^{-1}$)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
ylabel({'Resources, $R^{*}$ ($\mu$g/ml)'},...
    'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',20'); 
title('Invasion fitness with $0 < p < 1$, $\mathcal{R}_{0}$',...
    'fontsize',20,'interpreter','latex');


figure(5);
%heatmap of horizontal fitness
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
tmpt1 = text(1e7, 1.9, {'$\mathcal{R}_{0}>1$'});
set(tmpt1,'fontsize',25,'interpreter','latex');
tmpt2 = text(2e5, 0.8, {'$\mathcal{R}_{0} < 1$'});
set(tmpt2,'fontsize',25,'interpreter','latex');
x = [0.30 0.26]; % start
y = [0.25 0.18]; % end
annotation('arrow',x,y, 'LineWidth',2, 'HeadLength',15, 'HeadWidth',15)

% lables, axises
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(round(min(R_set),1):0.2:round(max(R_set),1)); 
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

ff1 = figure(1);
ff1.Units = 'inches';
Width = 17; Height = 17;
ff1.PaperSize = [Width, Height];
ff1.PaperPosition = [0 0 Width, Height];
ff1.Position = [0 0 Width, Height];
Name1 = 'present_model_resource_explicit_horizontal_fitness';
print(ff1, Name1, '-depsc2','-r600');
print(ff1, Name1, '-dpdf','-r600');


% Output files
ff2 = figure(2);
ff2.Units = 'inches';
Width = 17; Height = 17;
ff2.PaperSize = [Width, Height];
ff2.PaperPosition = [0 0 Width, Height];
ff2.Position = [0 0 Width, Height];
Name2 = 'present_model_resource_explicit_vertical_fitness_direct_benefits';
print(ff2, Name2, '-depsc2','-r600');
print(ff2, Name2, '-dpdf','-r600');

% Output files
ff3 = figure(3);
ff3.Units = 'inches';
Width = 17; Height = 17;
ff3.PaperSize = [Width, Height];
ff3.PaperPosition = [0 0 Width, Height];
ff3.Position = [0 0 Width, Height];
Name3 = 'present_model_resource_explicit_lysogeny_indirect_benefits';
print(ff3, Name3, '-depsc2','-r600');
print(ff3, Name3, '-dpdf','-r600');

% Output files
ff4 = figure(4);
ff4.Units = 'inches';
Width = 17; Height = 17;
ff4.PaperSize = [Width, Height];
ff4.PaperPosition = [0 0 Width, Height];
ff4.Position = [0 0 Width, Height];
Name4 = 'present_model_resource_explicit_intermed_fitness_indirect_benefits';
print(ff4, Name4, '-depsc2','-r600');
print(ff4, Name4, '-dpdf','-r600');

ff5 = figure(5);
ff5.Units = 'inches';
Width = 17; Height = 17;
ff5.PaperSize = [Width, Height];
ff5.PaperPosition = [0 0 Width, Height];
ff5.Position = [0 0 Width, Height];
Name5 = 'present_model_resource_explicit_intermed_fitness_direct_benefits';
print(ff5, Name5, '-depsc2','-r600');
print(ff5, Name5, '-dpdf','-r600');

