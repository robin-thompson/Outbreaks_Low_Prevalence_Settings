function [ext_sim,ext_kol,q,t,R0_inv,t_start_vec] = outbreak_func(t_final1,run_sim,run_prob,run_q,Thresh,sim,t_run,beta,mu,tv,inc,tR)

%%% Initialise variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_start_vec = [0:inc:t_final1];
lx = length(t_start_vec);
R0_inv = mu./beta;
ll = length(R0_inv);

%% Calculate SOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ext_sim = ones(lx,1);
init = 1;        % Number of initial infectives
if run_sim == 1
for xx = 1:lx
co = 0;
t_start = t_start_vec(xx);

for k=1:sim  % Run in parallel
tt = 0;
I = 0;
R = 0;
tt(1)=t_start;
I(1)=init;
R(1) = 0;
j=1;
S = 1;
while I(j)>0 && S <= Thresh % Run simulation until extinction or if threshold is reached
         u1 = rand; % uniform random number
         u2 = rand; % uniform random number
         if tt(j) < tv(length(tv))
            beta_t = interp1(tv,beta,tt(j));
         else
            beta_t = beta(length(beta));
         end 
         % Calculate time jump and probability of infection or recovery
         if tt(j) < tv(length(tv)) || beta_t == beta(length(beta))
         a = beta_t*I(j)+(mu)*I(j);
         probi = beta_t/(beta_t+mu);
         tt(j+1)=tt(j)-log(u1)/a;
         else
         a_fun = integral(@(s) (beta_t(s)+ mu)*I(j),tt(j),tt(j)+tau)
         u = fminsearch(a_fun(tau)+log(u1),1)
         probi = beta_t/(beta_t+mu);
         tt(j+1)=tt(j) + u;
         end
    if u2 <= probi
            I(j+1)=I(j)+1;  % infection event
            R(j+1)=R(j);
    else
            I(j+1)=I(j)-1;  % recovery event        
            R(j+1)=R(j)+1;
    end
j=j+1; 
S = I(j);

end

if I(length(I)) ~= 0    % Determine if threshold is reached and if an outbreak occured
co = co+1;
end

end
ext_sim(xx) = 1 - co/sim;   % Calculate the proportion of simulations that result in an outbreak

end

end

%% Call NOR Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ext_kol = ones(lx,1);
if run_prob == 1
for xx = 1:lx  % Run in parallel 
t_start = t_start_vec(xx);
ext_kol(xx) =  NOR_func(t_start,Thresh,t_run,beta,tv,mu);
end
end

%% Calculate COR & IOR %%%%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(tR)
    sett = 4;
    q_end = 0;
else
    sett = 2;
    q_end = 25;
end


t = [];
q = [];

if run_q == 1
options = odeset('RelTol',1e-15,'AbsTol',1e-16);

if sett == 1 || sett == 2
if q_end == 0

if sett == 1    
y0 = ext_sim(1);
elseif sett == 2
y0 = NOR_func(q_end,100,t_run,beta,tv,mu);
end

tspan2 = linspace(q_end,t_final1);
[t,q] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan2,y0,options);  %% Use ODE15s for accuracy as stiff problem  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif q_end == t_final1
if sett == 1    
y0 = ext_sim(q_end + 1);
elseif sett == 2
end

tspan2 = linspace(t_final1,0);
[t,q] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan2,y0,options);  %% Use ODE15s for accuracy as stiff problem  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else    
% q equation: solve backwards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan1 = linspace(q_end,0);

if sett == 1    
y0 = ext_sim(q_end + 1);
elseif sett == 2
y0 = NOR_func(q_end,Thresh,t_run,beta,tv,mu);
end
[t1,rho1] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan1,y0,options);  %% Use ODE15s for accuracy as stiff problem

% q equation: solve forwards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = rho1(1);
tspan2 = linspace(q_end,t_final1);
[t2,rho2] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan2,y0,options);  %% Use ODE15s for accuracy as stiff problem
t1(1) = [];
rho1(1) = [];
t = [flip(t1);t2];
q = [flip(rho1);rho2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif sett == 3
% q equation: solve backwards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan1 = linspace(tR,0);
y0 = 1 - 1/Thresh;                 
[t,q] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan1,y0,options);  %% Use ODE15s for accuracy as stiff problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif sett == 4
% q equation: solve backwards %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan1 = linspace(t_final1,0);
y0 = R0_inv(ll);
[t,q] = ode15s(@(t,Y) COR_ode(t,Y,beta,tv,mu),tspan1,y0,options);  %% Use ODE15s for accuracy as stiff problem
end

end

end





