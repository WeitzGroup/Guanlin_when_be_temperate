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
alpha_s_set = -1:1e-2:(1 - p.rho_SL/p.b_SL) - 4e-1; 
% theta_set = 1e-2:1e-2:10; % auxillary free var, theta = C - R*(alpha_s)
theta_set = 2e-2:1e-2:50;
R_set_SL = p.Rmod_SL.*p.rho_SL./((1 - alpha_s_set).*p.b_SL - p.rho_SL);
S_set_SL = theta_set./p.e_SL;
Rhor_set = zeros(length(S_set_SL),1);
% Rhor vs. Rver on (R*, S*) plane

data_max_strategy = zeros(length(S_set_SL),length(R_set_SL));
% data = maximal reproduction numbers
data_max_fitness = zeros(length(S_set_SL),length(R_set_SL));
C_matrix = zeros(length(S_set_SL),length(R_set_SL));

for i = 1:length(S_set_SL)
    for j = 1:length(R_set_SL)
        S_SL = S_set_SL(i); R_SL = R_set_SL(j);
        C_matrix(i,j) = theta_set(i) + R_SL;
        psi_mod = p.b_SL*R_SL/(R_SL + p.Rmod_SL);
        % two exclusive strategy fitness
        Rhor_SL = p.phi_SL*S_SL*p.beta_SL/p.rho_SL;
        Rver_SL = psi_mod/(p.rho_SL + p.tau_SL);
        data_max_fitness(i,j) = max(Rhor_SL, ((p.rho_SL + p.tau_SL)/(p.rho_SL + p.tau_SL + p.gamma_min))*Rver_SL);
        Rhor_set(i) = Rhor_SL;
        if max(Rhor_SL, ((p.rho_SL + p.tau_SL)/(p.rho_SL + p.tau_SL + p.gamma_min))*Rver_SL) > 1
            feasible_strategy_flag = true;
        else
            feasible_strategy_flag = false;
        end
        
        if feasible_strategy_flag % if resident strategy can invade
          
           % classify different cases in robustness
           if Rhor_SL >= Rver_SL
               data_max_strategy(i,j) = 1;
           else
               data_max_strategy(i,j) = 2;              
           end       
       else % Not feasible
           data_max_strategy(i,j) = 0;
        end
       
    end
end

S_low = S_set_SL(5); S_high = S_set_SL(end - 2000); 
R_low = R_set_SL(20); R_high = R_set_SL(end - 5); 
envir_temp = [R_high, S_low]; envir_lytic = [R_low, S_high]; 
R_c = p.Rmod_SL./(p.b_SL./(Rhor_set.*(p.rho_SL + p.tau_SL)) - 1); % transition curve
figure(1);
maph = [ 0.8 0.8 0.8 % gray
         1 0 0 % red
         0 0 1]; % blue   
% plot strategy transition on (S*, R*) plane 
h = pcolor(S_set_SL,R_set_SL,data_max_strategy'); hold on;
set(h, 'EdgeColor', 'none'); hold on;
colormap(maph); alpha(h,0.7);
C=caxis; caxis([C(1),C(2)]); % lock down the current colormap

levels =  [0.8, 1.2];
[M,c] = contour(S_set_SL,R_set_SL, data_max_fitness', levels, 'k','ShowText', 'on');
c.LineWidth = 1; 
clabel(M,c,'FontSize',18);
clabel(M,c,'LabelSpacing',300);
hold on;

levels =  [2, 4, 8, 16, 32];
[M,c] = contour(S_set_SL,R_set_SL, data_max_fitness', levels, 'k','ShowText', 'on');
c.LineWidth = 1; 
clabel(M,c,'FontSize',18);
clabel(M,c,'LabelSpacing',100);
hold on;


plot(S_set_SL,R_c, '.k', 'MarkerSize',15); hold on;

%{
plot(S_low,R_high,...
    'dk','MarkerSize',22,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
plot(S_high,R_low, ...
    'sk','MarkerSize',22,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
%}
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20);
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); yticks(0.7:0.3:R_set_SL(end));
xlabel({'Susceptible population, $S^{*} (ml^{-1})$ '}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ylabel({'Resources, $R^{*} (\mu g/ml)$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 

figure(2); % R0(p,gamma | S_low, R_high), i.e., fitness vs. strategy 
colormap jet; 
p_set = 0:1e-2:1; gamma_set = 1e-2:1e-2:1;

subplot(1,2,1);
R0_SL = R_fun_strategy(p, envir_temp, p_set, gamma_set, 3);
s = surf(p_set, gamma_set, R0_SL','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(end),gamma_set(1),max(max(R0_SL)),'.k','MarkerSize',50);
% R0 = 1;
hold on;
plot3(ones(length(gamma_set),1).*p_set(end),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2);
hold on;
plot3(ones(length(gamma_set),1).*p_set(1),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2);
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(1),...,
    ones(length(p_set),1),'--k','LineWidth',2);
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(end),...,
    ones(length(p_set),1),'--k','LineWidth',2);
set(gca,'YScale','log'); set(gca,'ZScale','log'); zlim([1e-2 1e2]);
axis square; grid off; box on;
ax = gca; ax.BoxStyle = 'full'; view(-53,40);
ax = gca; ax.BoxStyle = 'full'; view(-53,40);
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;ax.ZAxis.LineWidth = 1.5; 
set(gca,'FontSize',20);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);



% R0(p,gamma | S_high, R_low), i.e., fitness vs. strategy 
figure(2);colormap jet; 
subplot(1,2,2);
R0_SL = R_fun_strategy(p, envir_lytic, p_set, gamma_set, 3);
s = surf(p_set, gamma_set, R0_SL','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(1).*ones(length(gamma_set),1),gamma_set,max(R0_SL),'k','LineWidth',10);
% R0 = 1;
hold on;
plot3(ones(length(gamma_set),1).*p_set(end),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2);
hold on;
plot3(ones(length(gamma_set),1).*p_set(1),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2);
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(1),...,
    ones(length(p_set),1),'--k','LineWidth',2);
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(end),...,
    ones(length(p_set),1),'--k','LineWidth',2);
set(gca,'YScale','log'); set(gca,'ZScale','log'); zlim([1e-2 1e2]);
axis square; grid off; box on;
ax = gca; ax.BoxStyle = 'full'; view(-53,40);
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;ax.ZAxis.LineWidth = 1.5; 
set(gca,'FontSize',20);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);

