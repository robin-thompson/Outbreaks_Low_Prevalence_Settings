function [ext,p,t] = NOR_func(t_start,Thresh,t_run,beta,tv,mu)

p0=zeros(Thresh+1,1);
p0(2,1)=1;  % Start with one individual initially infected
tspan = linspace(t_start,t_start+t_run);

[t,p] = ode45(@(t,Y) my_ode(t,Y,beta,tv,mu,Thresh),tspan,p0);  % Solve forward kolmogorov equations
ext = p(size(p,1),1);                                          % Determine if extinction occured 

function dYdt = my_ode(t,Y,beta,tv,mu,Thresh)

% Forward kolmogorov equations
dYdt = T_func(t,beta,mu,tv,Thresh)*Y;

%%%% Build generator matrix Q
function T = T_func(t,beta,mu,tv,N)
T=zeros(N+1,N+1);
v=linspace(0,N,N+1);
if t < tv(length(tv))
    bt = interp1(tv,beta,t)*v;
else
    bt = beta(length(beta))*v;
end  
dt=mu*v;
dt(N+1) = 0;    % remove recoveries for N

  for i=2:N % Define the transition matrix
      T(i,i)=-bt(i)-dt(i); % diagonal entries
      T(i,i+1)=dt(i+1); % superdiagonal entries
      T(i+1,i)=bt(i); % subdiagonal entries
  end
T(1,2)=dt(2);
T(1,1) = 0;  
end
end

end