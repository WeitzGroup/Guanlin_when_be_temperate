p_set = 0:0.02:1; gamma_set = 1e-2:1e-2:1; % strategy space % Fixed parameters
p.r1 = 1.2;  % max growth rate of lysogenic cells
p.d = 0.2; % death rate of susceptible cells
p.d1 = 0.2; % death rate of infected cells
p.d2 = 0.2; % death rate of lysogenic cells
p.d3 = 0.2; % death rate of lytic cells
p.m = 0.4; % decay rate of phage
p.K = 2e8; % carrying capacity 
p.phi = 3.4e-10;  % adsorption rate
p.alpha = 2; % transit decision rate
p.eta = 1; % lysis rate
p.beta = 50; % burst rate

r_set = [1+5e-4:2e-4:1+1e-2, 1+2e-2:1e-2:1+1e-1, 1+2e-1:1e-1:2, 3:1:100].*p.d;
S_set = p.K.*(1 - p.d./r_set);
S_low = S_set(40); S_high = S_set(end-95); 

% Rhor and Rver 
Rhor = (p.beta*p.eta*p.phi*p.alpha).*S_set./...
    ((p.eta + p.d3).*(p.phi.*S_set + p.m).*(p.alpha + p.d1));
Rver = (p.r1/p.d2).*(1 - S_set./p.K);
loglog(S_set,Rhor,'r','LineWidth',4); hold on;
loglog(S_set,Rver,'b','LineWidth',4); hold on;
temp = 1e-2:1e-2:1e2; S_c = S_set(min(find(Rhor > Rver)) - 1);
loglog(S_c.*ones(length(temp),1), temp, 'k','LineWidth',3);
loglog(S_set,ones(length(S_set),1),'k--','LineWidth',2); hold on;
%{
plot(S_low,1e-2,...
    'dk','MarkerSize',20,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
plot(S_high,1e-2, ...
    'sk','MarkerSize',20,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','k'); hold on;
%}
axis square; 
ax = gca; 
xlim([min(S_set) max(S_set)]); ylim([1e-2 1e2]); 
set(gca,'FontSize',20);
xlabel({'Susceptible population, $S^{*} (ml^{-1})$ '}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ylabel({'Reproduction number'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;
xticks([1e5 1e6 1e7 1e8]); yticks([1e-2 1e-1 1e0 1e1 1e2]);
legend({'$\mathcal{R}^{hor}$','$\mathcal{R}^{ver}$'},'Interpreter','latex','FontSize',25); legend boxoff; 


figure(2); 
subplot(1,2,1);colormap jet;
R0 = R_fun_strategy(p, S_low, p_set, gamma_set, 1);
s = surf(p_set, gamma_set, R0','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(end),gamma_set(1),max(max(R0)),'.k','MarkerSize',50);
hold on;
% R0 = 1;
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
set(gca,'YScale','log');  set(gca,'ZScale','log'); zlim([1e-2 1e2]);
axis square; grid off; box on;
ax = gca; ax.BoxStyle = 'full'; 
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;ax.ZAxis.LineWidth = 1.5; 
set(gca,'FontSize',20);
view(-53,40);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);

% R0(p,gamma | S_high), i.e., fitness vs. strategy 
% figure(1);
figure(2); 
subplot(1,2,2); colormap jet;
R0 = R_fun_strategy(p, S_high, p_set, gamma_set, 1);
s = surf(p_set, gamma_set, R0','FaceAlpha',0.5); s.EdgeColor = 'none';
hold on;
plot3(p_set(1).*ones(length(gamma_set),1),gamma_set,max(R0),'k','LineWidth',10);
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
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5; ax.ZAxis.LineWidth = 1.5; 
set(gca,'FontSize',20);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);


