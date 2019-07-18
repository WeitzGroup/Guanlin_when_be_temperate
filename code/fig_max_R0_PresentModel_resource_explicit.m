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
p.m = 0.4;

% assume gamma from [1e-2, 1];
p.gamma_min = 1e-2; p.gamma_max = 1;
theta_set = [0.01:0.01:0.1, 0.2:0.05:1, 2:0.5:15]; % augumented parameter range
alphaS_set = -1:0.02:0.4; % alphaS < 0, S grows faster than L; alphaS > 0, ...

% (R*, S*) is varied by the parameter space (alphaS, J);
data_alphaS = zeros(length(alphaS_set),length(theta_set));
data_J = zeros(length(alphaS_set),length(theta_set));

% Range of Susceptible hosts and resources
S_set = theta_set./(p.e*p.ds);
R_set = p.Rmod*p.ds./((1 - alphaS_set).*p.b0 - p.ds);

% Using S* to compute the range of horizonral fitness
R_hor_set = p.beta*p.eta*p.phi.*S_set.*p.alpha./...
             ((p.eta + p.di)*(p.phi*S_set + p.m)*(p.alpha + p.de));

% Transition curve, i.e., R_hor(S*) = R_ver(R*)
% R_c = transfun(S)
R_c = p.Rmod./(p.b0./(R_hor_set.*p.dl) - 1);

% Save horizontal and vertical transmission fitness at virus-free
data_Rver = zeros(length(alphaS_set),length(theta_set));
data_Rhor = zeros(length(alphaS_set),length(theta_set));

% data_max_strategy: 0 - not feasible, 1 - lytic (resident ),... 
% 2 - lysogenic(resident)
data_max_strategy = zeros(length(alphaS_set),length(theta_set));
% data = maximal reproduction numbers
data_max_fitness = zeros(length(alphaS_set),length(theta_set));

%% Task 1: Partition of R-S plane
%  Proceduere:
%  1. For each input parameter pair, (alphaS,J), compute two modes of fitness values
%  2. Identify the invasion strategy that maximizes R0 using step 1.
%  3. compute the R0 of maximal strategy if
%    R0 > 1, otherwise, set it gray (no strategy can invade given this environment)

feasible_strategy_flag = true; % initialize the feasible strategy, can invade virus-free

for i = 1:length(alphaS_set)
    for j = 1:length(theta_set)
        % assign two model parameters
        p.alphaS = alphaS_set(i);
        p.J = theta_set(j) + p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        % Save the parameter pair
        data_alphaS(i,j) = p.alphaS; data_J(i,j) = p.J;
        
        % compute virus-free environmental states
        R_star = p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        S_star = theta_set(j)/(p.e*p.ds);

        % compute two transmission fitness values
        R_hor = p.beta*p.eta*p.phi*S_star*p.alpha/...
             ((p.eta + p.di)*(p.phi*S_star + p.m)*(p.alpha + p.de));
        psi_star = p.b0*R_star/(R_star + p.Rmod);         
        R_ver = psi_star/p.dl;
        % save two fitness values
        data_Rver(i,j) = R_ver; data_Rhor(i,j) = R_hor;
        data_max_fitness(i,j) = max(R_hor, (p.dl/(p.dl + p.gamma_min))*R_ver);
        
        if max(R_hor, (p.dl/(p.dl + p.gamma_min))*R_ver) > 1
            feasible_strategy_flag = true;
        else
            feasible_strategy_flag = false;
        end
        
       if feasible_strategy_flag % if resident strategy can invade
          
           % classify different cases in robustness
           if R_hor > R_ver
               data_max_strategy(i,j) = 1;
           else
               data_max_strategy(i,j) = 2;              
           end       
       else % Not feasible
           data_max_strategy(i,j) = 0;
       end
       
    end
    
end

S_low = S_set(5); S_high = S_set(end - 10); 
R_low = R_set(10); R_high = R_set(end - 5); 
envir_temp = [R_high, S_low]; envir_lytic = [R_low, S_high]; 
R_c = p.Rmod./(p.b0./(R_hor_set.*p.dl) - 1); % transition curve
maph = [ 0.8 0.8 0.8 % gray
         1 0 0 % red
         0 0 1]; % blue   
     
figure(1);
pcolor(S_set,R_set,data_max_strategy); 
h = pcolor(S_set,R_set,data_max_strategy); hold on;
set(h, 'EdgeColor', 'none'); colormap(maph); hold on; alpha(h,0.7);
set(gca,'XScale','log'); axis square;
set(gca,'FontSize',20);
ax = gca; ax.XAxis.LineWidth = 2; ax.YAxis.LineWidth = 2;
xticks([1e5 1e6 1e7 1e8]); 
yticks(0.5:0.2:max(R_set)); 

C=caxis; caxis([C(1),C(2)]); % lock down the current colormap

levels =  [0.6, 0.8];
[M,c] = contour(S_set,R_set, data_max_fitness, levels, 'k','ShowText', 'on');
c.LineWidth = 1; 
clabel(M,c,'FontSize',18);
clabel(M,c,'LabelSpacing',80);
hold on;

levels =  [1.1, 1.5, 2, 3, 4];
[M,c] = contour(S_set,R_set, data_max_fitness, levels, 'k','ShowText', 'on');
c.LineWidth = 1; 
clabel(M,c,'FontSize',18);
clabel(M,c,'LabelSpacing',180);
hold on;

plot(S_set, R_c, 'k', 'LineWidth',4); hold on;
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
xlabel({'Susceptible population, $S^{*} (ml^{-1})$ '}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 
ylabel({'Resources, $R^{*} (\mu g/ml)$'}, 'Interpreter','latex', 'FontName', 'Times New Roman','FontSize',25'); 


%% Task 2: R0 vs. (p, gamma) - purely lysogenic favored
p_set = 0:1e-2:1;
gamma_set = p.gamma_min:1e-2:p.gamma_max;

figure(2);
subplot(1,2,1);
R0 = R_fun_strategy(p, envir_temp, p_set, gamma_set, 4);
s = surf(p_set, gamma_set, R0','FaceAlpha',0.5); s.EdgeColor = 'none';
colormap jet; 
hold on;
plot3(p_set(end),gamma_set(1),max(max(R0)),'.k','MarkerSize',50); % R0 = 1;
hold on;
plot3(ones(length(gamma_set),1).*p_set(end),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2); % R0 = 1;
hold on;
plot3(ones(length(gamma_set),1).*p_set(1),gamma_set(1:end),...,
    ones(length(gamma_set),1),'--k','LineWidth',2); % R0 = 1;
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(1),...,
    ones(length(p_set),1),'--k','LineWidth',2); % R0 = 1;
hold on;
plot3(p_set(1:end),ones(length(p_set),1).*gamma_set(end),...,
    ones(length(p_set),1),'--k','LineWidth',2); % R0 = 1;

set(gca,'YScale','log'); set(gca,'ZScale','log'); zlim([1e-2 1e2]);
axis square; grid off; box on;
ax = gca; ax.BoxStyle = 'full'; view(-53,40);
ax.XAxis.LineWidth = 1.5; ax.YAxis.LineWidth = 1.5;ax.ZAxis.LineWidth = 1.5; 
set(gca,'FontSize',20);
xlabel('$p$','Interpreter','latex','FontSize',30); 
ylabel('$\gamma$','Interpreter','latex','FontSize',30);
zlabel('$\mathcal{R}^{0}$','Interpreter','latex','FontSize',30);

figure(2);
subplot(1,2,2); 
R0 = R_fun_strategy(p, envir_lytic, p_set, gamma_set, 4);
s = surf(p_set, gamma_set, R0','FaceAlpha',0.5); s.EdgeColor = 'none';
colormap jet; 
hold on;
plot3(p_set(1).*ones(length(gamma_set),1),gamma_set,max(R0),'k','LineWidth',10);
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
