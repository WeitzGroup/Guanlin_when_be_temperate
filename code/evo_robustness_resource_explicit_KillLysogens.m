% Evolutionary invasion analysis on resource-explicit - SinkCells

%% Model parameters
% Note: Parameter space must lead to the *stable* endemic equilibirum

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

%alphaS_upp = 1 - p.ds/p.b0;
% assume gamma from [1e-2, 1];
p.gamma_min = 1e-2; p.gamma_max = 1;

%theta_set = [0.01, 0.1, 1, 10, 15];
%alphaS_set = [-1, -0.5, 0, 0.2, 0.4];

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

% robustness matrix:
% 0 - not feasible
% 1 - lytic (resident) robust
% -1 - lytic (resident) not robust, i.e., invaded by lysogenic
% 2 - lysogenic (resident) robust
% -2 - lysogenic (resident) not robust, i.e., invaded by lytic
data_evo_robustness = zeros(length(alphaS_set),length(theta_set));


%% ----------- Loop over the virus-free environmental states -----------
% Proceduere:
% 1. For each input parameter pair, (alphaS,J), compute two modes of fitness values
% 2. Identify the invasion strategy that maximizes R0 using step 1.
% 3. compute the R0 of maximal strategy and set it as resident strategy if
%    R0 > 1, otherwise, set it gray (no strategy can invade given this environment)
% 4. Run the resident dynamics and set the other extreme strategy as the
%    mutant strategy, compute invasion fitness of mutant strategy

feasible_strategy_flag = true; % initialize the feasible strategy, can invade virus-free

% numerical solver set-up
t0 = 0; tf = 2e3; tol = 1e-5; 
options = odeset('RelTol',tol,'AbsTol',tol);

for i = 1:length(alphaS_set)
    for j = 1:length(theta_set)
        disp([i,j]);
        % assign two model parameters
        p.alphaS = alphaS_set(i);
        p.J = theta_set(j) + p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        % Save the parameter pair
        data_alphaS(i,j) = p.alphaS; data_J(i,j) = p.J;
     
        % compute virus-free environmental states
        R_star = p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        S_star = theta_set(j)/(p.e*p.ds);
        % initial condition of system
        IC = [R_star,S_star,0,0,0,1e2,0,0,0,0]; % initial condition

        % compute two transmission fitness values
        R_hor = p.beta*p.eta*p.phi*S_star*p.alpha/...
             ((p.eta + p.di)*(p.phi*S_star + p.m)*(p.alpha + p.de));
        psi_star = p.b0*R_star/(R_star + p.Rmod);         
        R_ver = psi_star/p.dl;
        % save two fitness values
        data_Rver(i,j) = R_ver; data_Rhor(i,j) = R_hor;
        
        % set resident and mutant strategy
        if R_hor > R_ver % resident strategy is purely lytic
            p.pb = 0; p.gamma = p.gamma_max;
            p.pb_m = 1; p.gamma_m = p.gamma_min;
            % check if the maximal invasion possible
            res_inv_fitness = R_hor;
            if res_inv_fitness > 1
                feasible_strategy_flag = true;
            else
                feasible_strategy_flag = false;
            end
            
        elseif  R_hor < R_ver % resident strategy is purely lysogenic strategy
            p.pb = 1; p.gamma = p.gamma_min;
            p.pb_m = 0; p.gamma_m = p.gamma_max;
            % check if the maximal invasion possible
            res_inv_fitness = (p.dl/(p.dl + p.gamma))*R_ver;
            if res_inv_fitness > 1
                feasible_strategy_flag = true;
            else
                feasible_strategy_flag = false;
            end
            
        else % equal, not likely to happen, break out of loop
            break;
        end
        
       if feasible_strategy_flag % if resident strategy can invade
           % resident dynamics
           [t,x] = ode45(@ode_resource_explicit_KillLysogens,[t0,tf],IC,options,p);
           
           % compute invasion fitness of mutant strategy
           R_mut = calculateR0_KillLysogens(p, x(end,1:5), [p.pb_m, p.gamma_m]);
           
           % classify different cases in robustness
           if R_mut > 1 && R_hor > R_ver
               data_evo_robustness(i,j) = -1;
           elseif R_mut < 1 && R_hor > R_ver
               data_evo_robustness(i,j) = 1;
           elseif R_mut > 1 && R_hor < R_ver
               data_evo_robustness(i,j) = -2;
           elseif R_mut < 1 && R_hor < R_ver
               data_evo_robustness(i,j) = 2;
           end       
       else % Not feasible
           data_max_strategy(i,j) = 0;
           data_evo_robustness(i,j) = 0;
       end
       
    end
    
end

% write matrix into a text.file
fid = fopen('data_robustness_evo_robustness_resource_explicit_KillLysogens.txt','wt');
for ii = 1:size(data_evo_robustness,1)
    fprintf(fid,'%g\t',data_evo_robustness(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);

