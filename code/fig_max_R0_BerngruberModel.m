p_set = 0:0.02:1; gamma_set = 1e-2:1e-2:1; % strategy space % Fixed parameters
% Fixed parameters from TS-1 of paper
p.r1_B = 1.4; % growth rate of infected cells
p.K_B = 1e9; % carrying capacity
p.delta_B = 1; % fidelity of vertical transmission
p.beta_B = 200; % burst size
p.m_B = 0.75; % dilution rate
p.b_B = 1e-2; % probability of fusion after adsorption
p.phi_B = 1e-8; % adsorption constant
% Susceptibel host S* is varying in maximum growth rate r
r_set_B = [1+1e-4:1e-4:1+1e-2, 1+2e-2:1e-2:1+1e-1, 1 + 1e-1:1e-1:2, 3:1:50].*p.m_B;
%r_set_B = [1+2e-2:1e-2:1+1e-1, 1 + 2e-1:1e-1:2, 3:1:100].*p.m_B;
S_set_B = p.K_B.*(1 - p.m_B./r_set_B);
% Rhor and Rver 
Rhor_B = (p.b_B*p.beta_B*p.phi_B).*S_set_B./(p.phi_B.*S_set_B + p.m_B);
Rver_B = (p.r1_B.*p.delta_B/p.m_B).*(1 - S_set_B./p.K_B);

S_low_B = S_set_B(10); S_high_B = S_set_B(end-45); 

figure(1); % plot Rhor(S*) vs. Rver(S*)
colormap jet;
loglog(S_set_B,Rhor_B,'r','LineWidth',4); hold on;
loglog(S_set_B,Rver_B,'b','LineWidth',4); hold on;
temp = 1e-2:1e-2:1e2; S_c = S_set_B(min(find(Rhor_B > Rver_B)));
loglog(S_c.*ones(length(temp),1), temp, 'k','LineWidth',3);
%{
plot(S_low_B ,1e-2,...
    'dk','MarkerSize',20,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
plot(S_high_B ,1e-2, ...
    'sk','MarkerSize',20,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
%}
loglog(S_set_B,ones(length(S_set_B),1),'k--','LineWidth',2);

axis square; 
xlim([min(S_set_B) max(S_set_B)]); ylim([1e-2 1e2]); 
ax = gca;set(gca,'FontSize',20);
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
xlabel({'Susceptible population, $S^{*} (ml^{-1})$ '}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ylabel({'Reproduction number'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$'},'Interpreter','latex','FontSize',25); legend boxoff; 
   

figure(2); % R0(p,gamma | S_low), i.e., fitness vs. strategy 
colormap jet;
subplot(1,2,1);

R0_B = R_fun_strategy(p, S_low_B, p_set, gamma_set, 2);
s = surf(p_set, gamma_set, R0_B','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(end),gamma_set(1),max(max(R0_B)),'.k','MarkerSize',50);
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

% R0(p,gamma | S_high), i.e., fitness vs. strategy 
figure(2);
colormap jet;
subplot(1,2,2);
R0_B = R_fun_strategy(p, S_high_B, p_set, gamma_set, 2);
s = surf(p_set, gamma_set, R0_B','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(1).*ones(length(gamma_set),1),gamma_set,max(R0_B),'k','LineWidth',10);
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
set(gca,'YScale','log'); set(gca,'ZScale','log'); 
axis square; grid off; box on;
ax = gca; ax.BoxStyle = 'full'; view(-53,40);
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;ax.ZAxis.LineWidth = 1.5; 
zlim([1e-2 1e2]);
set(gca,'FontSize',20);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);