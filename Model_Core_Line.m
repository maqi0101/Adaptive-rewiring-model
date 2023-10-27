% This  code is the core program for adaptive rewiring model, 
% which is used to generate data and figures for line diagrams

clc;
clear; 

Sp=30; % species number of plant
Sh=30; % species number of herbivore
Sm=30; % species number of pollinator
mu_r=1; % used to regulate the intrinsic growth rate of species
r=mu_r*ones(Sp+Sh+Sm,1); % intrinsic growth rate of species
h=0.1; % Half-saturation constant for trophic and mutualistic interactions
e=0.8; % conversation coefficient of predation
%%
cM = 0.15;   % connectance of mutualistic or antagonistic subnetworks
trSpan = 1;  % length of niche axis [0,1]
nw=0.1;    % niche width
int_comp=0.01;    % strength of competition
% int_mut=0.1;   % strength of mutualism 
% int_ant=0.1;   % strength of antagonism
eta=1;   % the parameters used to adjust the probability of disconnection
%%
int_m=0.05:0.05:0.5;     % range of  mutualistic strength 
int_a=[0.05,0.15,0.25];  % range of  antagonistic strength 
%%
rep=60; % repetitions
ant_N=zeros(length(int_a),length(int_m),rep); % antagonistic nestedness 
mut_N=zeros(length(int_a),length(int_m),rep); % mutualistic nestedness
ant_Q=zeros(length(int_a),length(int_m),rep); % antagonistic modularity
mut_Q=zeros(length(int_a),length(int_m),rep); % mutualistic modularity
lam=zeros(length(int_a),length(int_m),rep);   % resilience of networks
DC_ph=cell(length(int_a),length(int_m),rep);  % the degree centrality of plants in antagonism
DC_pm=cell(length(int_a),length(int_m),rep);  % the degree centrality of plants in mutualism 
P_mut=cell(length(int_a),length(int_m),rep);   % per capita energy intake from mutualistic interactions
P_ant=cell(length(int_a),length(int_m),rep);   % per capita energy loss from antagonistic interactions
P_abudance=cell(length(int_a),length(int_m),rep);  % species biomass of plant
%%
m_tep=10000; % last 10^4 rewiring attempts.
%%
for kk=1:rep
    trait=rand(Sp+Sh+Sm,1); % the central position of each species
    beta_P = zeros(Sp);  % per capita competitive rate of plant 
    beta_H = zeros(Sh);  % per capita competitive rate of herbivore
    beta_M = zeros(Sm);  % per capita competitive rate of pollinator 
    %
    for i=1:Sp
        for j=i+1:Sp
            beta_P(i,j)=int_comp*overlap(trait(i),trait(j),nw,trSpan); % interspecific competition
            beta_P(j,i)=beta_P(i,j);
        end
    end
    beta_P = beta_P + eye(Sp); % intraspecific competition
    %
    for i=1:Sh
        for j=i+1:Sh
            beta_H(i,j)=int_comp*overlap(trait(Sp+i),trait(Sp+j),nw,trSpan); % interspecific competition
            beta_H(j,i)=beta_H(i,j);
        end
    end
    beta_H=beta_H + eye(Sh); % intraspecific competition
    %
    for i=1:Sm
        for j=i+1:Sm
            beta_M(i,j)=int_comp*overlap(trait(Sp+Sh+i),trait(Sp+Sh+j),nw,trSpan); % interspecific competition
            beta_M(j,i)=beta_M(i,j);
        end
    end
    beta_M=beta_M + eye(Sm);% intraspecific competition
    %%%%%%%%%%%%%%%%%%%%%%
    alpha=zeros(Sp,Sh); % record the niche overlap of antagonistic interactions
    for i=1:Sp
        for j=1:Sh
            if rand < cM
                alpha(i,j)=overlap(trait(i),trait(Sp+j),nw,trSpan);
            end
        end
    end
    %
    gamma=zeros(Sp,Sm); % record the niche overlap of mutualistic interactions
    for i=1:Sp
        for j=1:Sm
            if rand<cM
                gamma(i,j)=overlap(trait(i),trait(Sp+Sh+j),nw,trSpan);  
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%  
    for ii=1:length(int_a)
        int_ant=int_a(ii); % strength of antagonism
        alpha_PH=int_ant*alpha; % per capita antagonistic rate 
        for jj=1:length(int_m)
            int_mut=int_m(jj);  % strength of mutualism
            gamma_PM=int_mut*gamma; % per capita mutualistic rate 
                    
            counter_max = 100000; % total number of timesteps
            counter=0;   % current count
            idstart = 0; % record the number of the animal
            Type = 0;    % rewiring in H or M guild     1:H->P; 2:M->P
            mu_N0 = 0.1; 
            N0 = mu_N0*ones(Sp+Sh+Sm,1);  % initial abundance
            to = 0;   % initial time
            tf = 50;  % integration time within each time step 
            
            % several topological properties for ecological networks of the last 10^4 rewiring attempts.
            ant_N_ijk=zeros(1,m_tep); % antagonistic nestedness
            ant_Q_ijk=zeros(1,m_tep); % antagonistic modularity
            mut_N_ijk=zeros(1,m_tep); % mutualistic nestedness
            mut_Q_ijk=zeros(1,m_tep); % mutualistic modularity
            lam_ijk=zeros(1,m_tep); % resilience
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % adaptive rewiring, updating niche proximities, integration
            while counter<counter_max
                counter = counter+1;
                
                % rewiring, updating niche proximities
                Type=randi(2);   % rewiring in H or M guild (1:H->P; 2:M->P)
                [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old_rw,int_old]=rewiring(Type,N0,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta);
                
                % integration
                [t y] = ode45(@(t,y)LV_PHM(t,y,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h),[to to+tf],N0);%%%重连后，运行一段时间
                to = t(end);
                N0=y(end,:)'; % set initial condition for next interval 
                
                % resilience
                if counter>counter_max-m_tep
                    jacob_mat=get_jacmat(N0,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h);
                    lam_ijk(counter+m_tep-counter_max)=-max(real(eig(jacob_mat)));
                end
                
                % link recovery if abundance does not increase
                if counter>1
                    if Type==1 
                        gain_new_rw=N0(idstart+Sp); % The biomass of this species after rewiring
                    else
                        gain_new_rw=N0(idstart+Sp+Sh);
                    end
                    if gain_old_rw > gain_new_rw     % wire back to previous partner
                        [alpha_PH,gamma_PM]=wire_back(alpha_PH,gamma_PM,int_old,idstart,idend_new,idend_old,Type);
                    end
                end
                
                % calculate nestedness and modularity
                if counter>counter_max-m_tep
                    % antagonistic nestedness and modularity
                    [nodf,qb,Nm] = cal_structure(alpha_PH);
                    ant_N_ijk(counter+m_tep-counter_max)=nodf;
                    ant_Q_ijk(counter+m_tep-counter_max)=qb;
                    % mutualisti nestedness and modularity
                    [nodf,qb,Nm] = cal_structure(gamma_PM);
                    mut_N_ijk(counter+m_tep-counter_max)=nodf;
                    mut_Q_ijk(counter+m_tep-counter_max)=qb;
                end
                fprintf('\nrepeat=%d, row = %d, list = %d,  count=%d\n\n', kk, ii, jj, counter);
            end
            
            P_abud=N0(1:Sp);   % abundance of plants in the final network
            H_abud=N0(Sp+1:Sp+Sh);  % abundance of herbivores in the final network
            M_abud=N0(Sp+Sh+1:Sp+Sh+Sm); % abundance of pollinator in the final network
            
            % interaction relationship; plant’ node degree 
            mat_pm=logical(gamma_PM);
            mat_ph=logical(alpha_PH);
            deg_pm_ijk=sum(mat_pm,2)';
            deg_ph_ijk=sum(mat_ph,2)';
            
            % energy flow situation
            mut_mat=(diag(M_abud)*gamma_PM')./(1+h*repmat((mat_pm*M_abud)',[Sm,1]));%check
            ant_mat=(alpha_PH*diag(H_abud))./(1+h*repmat((mat_ph'*P_abud)',[Sp,1]));%check
            
            % the average of the last 10^4 times
            ant_N(ii,jj,kk)=mean(ant_N_ijk);
            mut_N(ii,jj,kk)=mean(mut_N_ijk);
            ant_Q(ii,jj,kk)=mean(ant_Q_ijk);
            mut_Q(ii,jj,kk)=mean(mut_Q_ijk);
            lam(ii,jj,kk)=mean(lam_ijk);
            
            % degree centrality and energy flow
            DC_pm{ii,jj,kk}=deg_pm_ijk/Sm; % the degree centrality of plants in antagonism
            DC_ph{ii,jj,kk}=deg_ph_ijk/Sh; % the degree centrality of plants in mutualism 
            P_mut{ii,jj,kk}=sum(mut_mat);  % per capita energy intake from mutualistic interactions
            P_ant{ii,jj,kk}=sum(ant_mat,2)'; % per capita energy loss from antagonistic interactions
            P_abudance{ii,jj,kk}=P_abud;  % species biomass of plant
            
        end
    end
end
%%
% Calculate correlation
corr_DC=zeros(length(int_a),length(int_m),rep);  % correlation between correlations d_mut and d_mut
corr_bio=zeros(length(int_a),length(int_m),rep); % correlation between correlations b_mut and b_mut
corr_d_b=zeros(length(int_a),length(int_m),rep); % correlation between correlations (d_mut - d_mut) and (b_mut - b_mut)

for ii=1:length(int_a)
    for jj=1:length(int_m)
        for kk=1:rep
            dc_mut=DC_pm{ii,jj,kk};
            dc_ant=DC_ph{ii,jj,kk};
            diff_dc=dc_mut-dc_ant;
            bio_mut=P_mut{ii,jj,kk};
            bio_ant=P_ant{ii,jj,kk};
            diff_bio=bio_mut-bio_ant;
            cor1=corrcoef(dc_mut,dc_ant);
            cor2=corrcoef(bio_mut,bio_ant);
            cor3=corrcoef(diff_dc,diff_bio);
            corr_DC(ii,jj,kk)=cor1(2);
            corr_bio(ii,jj,kk)=cor2(2);
            corr_d_b(ii,jj,kk)=cor3(2);
        end
    end
end

%%
% the mean of repetitions
N_ant=mean(ant_N,3);
N_mut=mean(mut_N,3);
Q_ant=mean(ant_Q,3);
Q_mut=mean(mut_Q,3);
lambda=mean(lam,3);
dc_mean=mean(corr_DC,3);
bio_mean=mean(corr_bio,3);
dc_bio_mean=mean(corr_d_b,3);

% the standard deviation of repetitions
N_ant_err=std(ant_N,1,3);
N_mut_err=std(mut_N,1,3);
Q_ant_err=std(ant_Q,1,3);
Q_mut_err=std(mut_Q,1,3);
lambda_err=std(lam,1,3);
dc_err=std(corr_DC,1,3);
bio_err=std(corr_bio,1,3);
dc_bio_err=std(corr_d_b,1,3);
%%

figure(1)
errorbar(int_m,Q_ant(1,:),Q_ant_err(1,:));
hold on
errorbar(int_m,Q_ant(2,:),Q_ant_err(2,:));
errorbar(int_m,Q_ant(3,:),Q_ant_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Q_a_n_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%
figure(2)
errorbar(int_m,Q_mut(1,:),Q_mut_err(1,:));
hold on
errorbar(int_m,Q_mut(2,:),Q_mut_err(2,:));
errorbar(int_m,Q_mut(3,:),Q_mut_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Q_m_u_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%
figure(3)
errorbar(int_m,N_ant(1,:),N_ant_err(1,:));
hold on
errorbar(int_m,N_ant(2,:),N_ant_err(2,:));
errorbar(int_m,N_ant(3,:),N_ant_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('N_a_n_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%
figure(4)
errorbar(int_m,N_mut(1,:),N_mut_err(1,:));
hold on
errorbar(int_m,N_mut(2,:),N_mut_err(2,:));
errorbar(int_m,N_mut(3,:),N_mut_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('N_m_u_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%
figure(5)
errorbar(int_m,lambda(1,:),lambda_err(1,:));
hold on
errorbar(int_m,lambda(2,:),lambda_err(2,:));
errorbar(int_m,lambda(3,:),lambda_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Resilience')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%
figure(6)
errorbar(int_m,dc_mean(1,:),dc_err(1,:));
hold on
errorbar(int_m,dc_mean(2,:),dc_err(2,:));
errorbar(int_m,dc_mean(3,:),dc_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Correlation between d_m_u_t and d_a_n_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%%
figure(7)
errorbar(int_m,bio_mean(1,:),bio_err(1,:));
hold on
errorbar(int_m,bio_mean(2,:),bio_err(2,:));
errorbar(int_m,bio_mean(3,:),bio_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Correlation between b_m_u_t and b_a_n_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
%%
figure(8)
errorbar(int_m,dc_bio_mean(1,:),dc_bio_err(1,:));
hold on
errorbar(int_m,dc_bio_mean(2,:),dc_bio_err(2,:));
errorbar(int_m,dc_bio_mean(3,:),dc_bio_err(3,:));
hold off
xlabel('\Omega_m')
ylabel('Correlation between  d_m_u_t-d_a_n_t  and  b_m_u_t-b_a_n_t')
legend({'\Omega_p=0.05','\Omega_p=0.15','\Omega_p=0.25'})
title('\Omega_c=0.01')
