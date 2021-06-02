%%% Function to simulation figure 2 %%%%%%%%%%%%%%%
function [] = fig_2_func(S,mu,del,n,eta1,eta2,X_vec,R0m_vec)

for j = 1:length(X_vec) % Run for loop to create figures for IoM & Israel with R0 = 3, R0 = 5
X = X_vec(j);
R0m = R0m_vec(j); 

%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if X == 1
load('Data_IoM.mat')    % Load IoM vaccination data
N1 = 84500;             % IoM population
dat_end = 98/7;         % Number of weeks' worth of real vaccination data
elseif X == 2
load('Data_Israel.mat') % Load Israel vaccination data 
N1 = 8772800;           % Israel population
dat_end = 124/7;        % Number of weeks' worth of real vaccination data
end

%% Initialise variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
country_label = {'Isle of Man','Israel'};
V = [0.674,1.036,1.644,1.96];    % 50%, 70%, 90%, 95% Confidence interval bands for reproduction number
R_vec = [R0m,R0m-V(1)*S,R0m+V(1)*S,R0m-V(2)*S,R0m+V(2)*S,R0m-V(3)*S,R0m+V(3)*S,R0m-V(4)*S,R0m+V(4)*S];
t = dat(:,1)./7;
t = t';
lt = length(t);
V0 = dat(:,2);
V1 = dat(:,3);
V2 = dat(:,4);
scen = 1;

%% Compute R(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(R_vec)
R = R_vec(i);  
beta0 = R*mu;   % beta for unvaccinated 
alpha = beta0*(V0 + eta1*V1 + eta2*V2);
R0_end = alpha/(mu*N1);
R0_end([lt-(del-1):lt]) = [];
R0_beg = R*ones(del,1);
R0(:,i) = [R0_beg;R0_end];
end

%% Compute start/end dates of simulations and shaded areas %%%%%%%%%%%
t = datetime(2020,12,18) + caldays(0:n);
dat_start = datetime(2020,12,18);
if X == 1     
dat_end = datetime(2020,12,18) + caldays(114);
elseif X == 2
dat_end = datetime(2020,12,18) + caldays(124);
end
t_end = datetime(2020,12,18) + caldays(n);

%% Plot Vaccination Proportion of Population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if j == 1 || j == 3 % Only plot once for each country

figure
fig_settings % Figure settings
fill([dat_end t_end t_end dat_end],[eps eps 1 1],'yellow','FaceAlpha',0.2, 'Edgecolor','none')
hold on
p1 = plot(t,V0/N1,'b','linewidth',1.5);
p2 = plot(t,V1/N1,'r','linewidth',1.5);
p3 = plot(t,V2/N1,'k','linewidth',1.5);
title(country_label(X),'fontsize',24);
xlim(datetime([2020 2021],[12 8],[18 20]))
xticks(datetime(2021,[1:8],1))
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'})
xlabel('Date (2021)','fontsize',24)
ylabel('Proportion of Population','fontsize',24)
legend([p1,p2,p3],{'Unvaccinated','One dose ($V_1(t)$)','Two doses ($V_2(t)$)'},'fontsize',18,'location','northeast')
box on
shg
end

%% Plot Rv(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
fill([dat_end t_end t_end dat_end],[eps eps 6 6],'yellow','FaceAlpha',0.2, 'Edgecolor','none')
a = [0.8,0.6,0.45,0.3];
patch_col = [0 0 1];
fill([t fliplr(t)], [R0(:,8)' fliplr(R0(:,9)')],[1 0 0],'Facecolor',patch_col,'FaceAlpha',a(4),'EdgeColor','none')
fill([t fliplr(t)], [R0(:,6)' fliplr(R0(:,7)')],[1 0 0],'Facecolor',patch_col,'FaceAlpha',a(3),'EdgeColor','none')
fill([t fliplr(t)], [R0(:,4)' fliplr(R0(:,5)')],[1 0 0],'Facecolor',patch_col,'FaceAlpha',a(2),'EdgeColor','none')
fill([t fliplr(t)], [R0(:,2)' fliplr(R0(:,3)')],[1 0 0],'Facecolor',patch_col,'FaceAlpha',a(1),'EdgeColor','none')
plot(t,R0(:,1),'r','linewidth',1.5);
xlabel('Date (2021)','fontsize',24)
ylabel({'Reproduction Number $(R_V(t))$';},'fontsize',24)
% title(country_label(X),'fontsize',24);
if X == 1
title(['Isle of Man, Rv(0)=',num2str(R0m)],'fontsize',24);
elseif X == 2
title(['Israel, Rv(0)=',num2str(R0m)],'fontsize',24);
end
xlim(datetime([2020 2021],[12 8],[18 20]))
xticks(datetime(2021,[1:8],1))
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'})
yline(1,'k--','linewidth',1.5)
ylim([0 6])
box on
shg

end

end


