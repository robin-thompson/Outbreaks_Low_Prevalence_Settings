%%% Function to simulation figure 3 %%%%%%%%%%%%%%%
function [] = fig_3_func(S,mu,del,n,eta1,eta2,X_vec,R0m_vec,t_final1,Thresh,P,sim)

%% Parmeters %%%%%%%%%%%%%%%%%%%%%%%%%
inc = 1;        % Time increment in weeks when simulating SOR and NOR
run_SOR = 1;     % 1: run SOR
run_NOR = 1;     % 1: run NOR
run_IOR_COR = 1; % 1: run IOR & COR
t_run = 100;     % Number of weeks into the future NOR simulation are run for

for j = 1:length(X_vec)
X = X_vec(j)
R0m = R0m_vec(j)

%%% Vector of initial R values %%%%%%%%%%%%%%%%%%%%%%%%
if P == 1
R_vec = R0m;
else
R_vec = linspace(R0m-1.96*S,R0m+1.96*S,P);
R_norm = normpdf(R_vec,R0m,S)';
end


%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if X == 1
load('Data_IoM.mat')    % Load IoM vaccination data
N1 = 84500;             % IoM population
elseif X == 2
load('Data_Israel.mat') % Load Israel vaccination data 
N1 = 8772800;           % Israel population
end

%% Initialise variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = dat(:,1)./7;
t = t';
lt = length(t);
tv = t;
V0 = dat(:,2);
V1 = dat(:,3);
V2 = dat(:,4);
tt = linspace(0,t_final1,1000);
clear SER CER IER PER lCER

%% Compute R(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(R_vec)
R = R_vec(i);    
beta0 = R*mu;  
alpha = beta0*(V0 + eta1*V1 + eta2*V2);
beta_end = alpha/N1;
beta_end([lt-(del-1):lt]) = [];
beta_beg = beta0*ones(del,1);
beta = [beta_beg;beta_end];
rsim = run_SOR;

%%% Determine if R(t) crosses 1 %%%%%%%%%%%%%%%%%%%%%%%%
A = beta(length(beta))/mu;
if A > 1
tR = NaN;
else
tR = 1;
end

%% Call function which calculates outbreak risks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ext_SOR,ext_NOR,q,t,R0_inv] = outbreak_func(t_final1,rsim,run_NOR,run_IOR_COR,Thresh,sim,t_run,beta,mu,tv,inc,tR);

% Format outbreak risks vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SOR(i,:) = 1-ext_SOR;
IOR(i,:) = max(1-R0_inv,0);
NOR(i,:) = 1-ext_NOR;
if isnan(tR)
tCOR = flip(t);
Y = flip(1-q);
COR(i,:) = interp1(tCOR,Y,tv);
else
tCOR = t;
Y = 1-q;
COR(i,:) = interp1(tCOR,Y,tv);
end


end

if P == 1
SORa = SOR;
IORa = IOR;
CORa = COR;
NORa = NOR;
else
SORa = sum(R_norm.*SOR)/sum(R_norm);
IORa = sum(R_norm.*IOR)/sum(R_norm);
CORa = sum(R_norm.*COR)/sum(R_norm);
NORa = sum(R_norm.*NOR)/sum(R_norm);
end

%% Compute start/end dates of simulations and shaded areas %%%%%%%%%%%
t_end = datetime(2020,12,18) + caldays(n);
caldays_end = [114,124];
dat_end = datetime(2020,12,18) + caldays(caldays_end(X));
date_start = datetime(2020,12,18);
t_start_vec = date_start + caldays(0:7:n);
tv = date_start + caldays(0:1:n);

%% Plot Figure Probabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
fig_settings()
hold on
fill([dat_end t_end t_end dat_end],[eps eps 1 1],'yellow','FaceAlpha',0.2, 'Edgecolor','none')
p1 = scatter(t_start_vec,SORa,60,[0,0.6,0],'filled');
p2 = plot(t_start_vec,NORa,'r','linewidth',1.5);
p3 = plot(tv,CORa,'b--','linewidth',1.5);
p4 = plot(tv,IORa,'k','linewidth',1.5);
xlabel('Date of Introduction (2021)','fontsize',24)
ylabel('Outbreak Risk','fontsize',24)
ylim([0 1])
if X == 1
title(['Isle of Man, Rv(0)=',num2str(R0m)],'fontsize',24);
elseif X == 2
title(['Israel, Rv(0)=',num2str(R0m)],'fontsize',24);
end
xlim(datetime([2020 2021],[12 8],[18 20]))
xticks(datetime(2021,[1:8],1))
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'})
legend([p1,p4,p3,p2],{'SOR','IOR','COR','NOR'},'fontsize',22,'location','northeast')
box on
shg

end


end



