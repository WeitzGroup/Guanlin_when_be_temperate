% Comparison of optimal strategies in maximizing R0 and r
% note : basic reproduction number (R0), growth rate (r)

%% ------------------ Parameters -------------------------------
p.d = 0.2; % death rate of susceptible cells
p.d1 = 0.2; % death rate of infected cells
p.d2 = 0.2; % death rate of lysogenic cells
p.d3 = 0.2; % death rate of lytic cells
p.phi = 1e-9;  % adsorption rate
p.alpha = 2; % transit decision rate (from exposed state to fate-determined state)
p.K = 2e8; % carrying capacity
p.eta = 1; % lysis rate
p.beta = 50; % burst rate
p.m = 0.4; % decay rate
p.r1 = 1.2; % max growth rate of lysogeny cells

% S = K(1 - d/r), varying 'r' to modulate the susceptible density
r_set = [1+5e-4:1e-4:1+1e-2, 1+5e-3:5e-3:1+1e-1, 1+5e-2:5e-2:2, 3:1:100].*p.d;

%% ------------------ Search for maximization ---------------------
% Initialize Savers
optimalr_p_saver = zeros(length(r_set),1); 
optimalR_p_saver = zeros(length(r_set),1); 
optimalr_gamma_saver = zeros(length(r_set),1); 
optimalR_gamma_saver = zeros(length(r_set),1); 
% Check the uniqueness of maximizers
uqr_saver = zeros(length(r_set),1); 
uqR_saver = zeros(length(r_set),1); 
 
for j = 1:length(r_set)
p.r = r_set(j); p.S_ss = p.K*(1 - p.d/p.r); % susceptible population
% initialize decision space
p_set = 0:0.02:1; gamma_set = 1e-3:5e-3:1e-1; 
r_saver = zeros(length(p_set),length(gamma_set)); % growth rate r(p,gamma)
R_saver = zeros(length(p_set),length(gamma_set)); % basic reproduction number R0(p,gamma)

% compute R_hor(S) and R_ver(S) 
%{ 
R_hor = p.beta*p.eta*p.phi*p.S_ss*p.alpha/...
    ((p.eta + p.d3)*(p.phi*p.S_ss + p.m)*(p.alpha + p.d1));
R_ver = (p.r1/p.d2)*(1 - p.S_ss/p.K);
%}
for k = 1:length(p_set)
    for i = 1:length(gamma_set)
    p.pb = p_set(k); p.gamma = gamma_set(i);
    % transmission matrix T:
    T = zeros(4); T(1,4) = p.phi*p.S_ss; T(2,2) = p.r1*(1 - p.S_ss/p.K);
    % transition matrix Sigma;
    Sigma = zeros(4); Sigma(1,1) = -(p.alpha + p.d1); Sigma(2,2) = -(p.gamma + p.d2); 
    Sigma(3,3) = -(p.eta + p.d3);  Sigma(4,4) = -(p.phi*p.S_ss + p.m); Sigma(2,1) = p.alpha*p.pb; 
    Sigma(3,1) = p.alpha*(1 - p.pb); Sigma(3,2) = p.gamma;  Sigma(4,3) = p.beta*p.eta;
    % compute r(p,gamma) and R0(p,gamma)
    r_saver(k,i) = max(real(eig(T + Sigma))); 
    R_saver(k,i) = max(eig(- T*(Sigma\eye(4)))); % saver
    end
end

[MR,IR] = max(R_saver(:)); [IR_row, IR_col] = ind2sub(size(R_saver),IR);
[Mr,Ir] = max(r_saver(:)); [Ir_row, Ir_col] = ind2sub(size(r_saver),Ir);
optimalr_p_saver(j) = p_set(Ir_row);optimalR_p_saver(j) = p_set(IR_row);
optimalr_gamma_saver(j) = gamma_set(Ir_col); optimalR_gamma_saver(j) = gamma_set(IR_col);

% uniqueness;
uqR_saver(j) = sum(sum(R_saver == max(max(R_saver))));
uqr_saver(j) = sum(sum(r_saver == max(max(r_saver))));
end

figure(1);
subplot(1,2,1);
yyaxis left;
semilogx(p.K.*(1 - p.d./r_set), optimalR_p_saver,'.b','MarkerSize',20); 
yticks([0 1]); yticklabels({'p^{*} = 0','p^{*} = 1'}); hold on;
yyaxis right;
semilogx(p.K*(1 - p.d./r_set((find(uqR_saver == 1)))),...
    optimalR_gamma_saver(find(uqR_saver == 1)),'.r','MarkerSize',20);
yticks([min(gamma_set) max(gamma_set)]); 
yticklabels({'\gamma^{*} = \gamma_{min}','\gamma^{*} = \gamma_{max}'});
xlabel('$\textbf{Susceptible population,}\ S^{*}$','Interpreter','latex','FontSize',15); 
yyaxis left; 
ylabel('$\textbf{lysogeny probability}$','Interpreter','latex','FontSize',15);
yyaxis right; 
ylabel('$\textbf{induction rate}$','Interpreter','latex','FontSize',15);
ylim([min(gamma_set) max(gamma_set)]); axis square; 
xlim([min(p.K.*(1 - p.d./r_set)) max(p.K.*(1 - p.d./r_set))]);
low_limit = min(optimalR_p_saver);
high_limit = max(optimalR_p_saver);
regime1 = p.K*(1 - p.d./r_set(1:min(find(optimalR_p_saver ~= 1))));
plotshaded(regime1,...
   [low_limit.*ones(1,length(regime1));...
   high_limit.*ones(1,length(regime1))],'g'); hold on;
regime2 = p.K*(1 - p.d./r_set(min(find(optimalR_p_saver ~= 1)):end));
plotshaded(regime2,...
   [low_limit.*ones(1,length(regime2));...
   high_limit.*ones(1,length(regime2))],'k'); hold on;

figure(1);
subplot(1,2,2);
yyaxis left;
semilogx(p.K.*(1 - p.d./r_set), optimalr_p_saver,'.b','MarkerSize',20); 
yticks([0 1]); yticklabels({'p^{*} = 0','p^{*} = 1'});
yyaxis right;
semilogx(p.K*(1 - p.d./r_set((find(uqr_saver == 1)))),...
    optimalr_gamma_saver((find(uqr_saver == 1))),'.r','MarkerSize',20);
yticks([min(gamma_set) max(gamma_set)]); 
yticklabels({'\gamma^{*} = \gamma_{min}','\gamma^{*} = \gamma_{max}'});
xlabel('$\textbf{Susceptible population,}\ S^{*}$','Interpreter','latex','FontSize',15); 
yyaxis left; 
ylabel('$\textbf{lysogeny probability}$','Interpreter','latex','FontSize',15);
yyaxis right; 
ylabel('$\textbf{induction rate}$','Interpreter','latex','FontSize',15);
ylim([min(gamma_set) max(gamma_set)]); axis square; 
xlim([min(p.K.*(1 - p.d./r_set)) max(p.K.*(1 - p.d./r_set))]);
low_limit = min(optimalr_p_saver);
high_limit = max(optimalr_p_saver);
regime1 = p.K*(1 - p.d./r_set(1:min(find(optimalr_p_saver ~= 1))));
plotshaded(regime1,...
    [low_limit.*ones(1,length(regime1));...
    high_limit.*ones(1,length(regime1))],'g'); hold on;
regime2 = p.K*(1 - p.d./r_set(min(find(optimalr_p_saver ~= 1)):end));
plotshaded(regime2,...
    [low_limit.*ones(1,length(regime2));...
    high_limit.*ones(1,length(regime2))],'k'); hold on;


