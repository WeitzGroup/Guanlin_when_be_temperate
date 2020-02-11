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

% fix induction rate 
p.gamma = 1e-2;

% (R*, S*) is varied by the parameter space (alphaS, J);
theta_set = [0.01:0.01:0.1, 0.2:0.05:1, 2:0.5:10]; % augumented parameter range
alphaS_set = -1:0.02:0.5; % alphaS < 0, S grows faster than L; alphaS > 0, ...

dl_set = [0.8, 0.08]; % indirect benefits and direct benefits
p.dl = dl_set(1); % indirect benefits

% saver matrix
equilibrium_id = zeros(length(alphaS_set),length(theta_set));
R_ver_endemic_saver = zeros(length(alphaS_set),length(theta_set));

for i = 1:length(alphaS_set)
    for j = 1:length(theta_set)
        % assign two model parameters
        p.alphaS = alphaS_set(i);
        p.J = theta_set(j) + p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        
        % compute virus-free environmental states
        R_star = p.Rmod*p.ds/((1 - p.alphaS)*p.b0 - p.ds);
        S_star = theta_set(j)/(p.e*p.ds);

        % compute two transmission fitness values
        R_hor = p.beta*p.eta*p.phi*S_star*p.alpha/...
             ((p.eta + p.di)*(p.phi*S_star + p.m)*(p.alpha + p.de));
        psi_star = p.b0*R_star/(R_star + p.Rmod);
        R_ver_ind = psi_star/p.dl; 
        
        if R_hor > 1
            % simulation 
            V_init = 1; t0 = 0; tf = 5e3;
            tol = 1e-6; options = odeset('RelTol',tol,'AbsTol',tol);
            IC = [R_star,S_star,0,0,0,V_init]; % initial condition
            p.pb = 0; p.gamma = 0.1; % gamma doesn't matter for lytic viruses
            [t,x] = ode45(@ode_resource_explicit_NoSuperinfection,[t0,tf],IC,options,p);
            if var(x(find(t > tf - 1e3), 1))/mean(x(find(t > tf - 1e3), 1)) < 1e-8
                equilibrium_id(i,j) = 2; % stable equilibrium
                endemic_R = x(end, 1); endemic_S = x(end, 2); % endemic environment at equilibrium
                psi_star_endemic = p.b0*endemic_R/(endemic_R + p.Rmod);      
                R_ver_ind_endemic = psi_star_endemic/p.dl; 
                R_ver_endemic_saver(i,j) = R_ver_ind_endemic;
            else % oscillatory equilibrium
                equilibrium_id(i,j) = 1; 
                R_ver_endemic_saver(i,j) = -1;
            end
        end
        
        disp([i,j]);
                      
    end
end

% write matrix into a text.file
fid = fopen('data_equilibrium_id_endemic.txt','wt');
for ii = 1:size(equilibrium_id,1)
    fprintf(fid,'%g\t',equilibrium_id(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);

% write matrix into a text.file
fid = fopen('data_vertical_fitness_endemic.txt','wt');
for ii = 1:size(R_ver_endemic_saver,1)
    fprintf(fid,'%g\t',R_ver_endemic_saver(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


