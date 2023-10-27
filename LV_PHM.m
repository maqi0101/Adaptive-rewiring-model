% This function is used to describe systems of differential equations

function dydt = LV_PHM(t,y,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h)
%%
dydt=zeros(Sp+Sh+Sm,1); % Sp,Sh,Sm:species number of plant, herbivore or pollinator
Pi=y(1:Sp); % plant 
Hi=y(Sp+1:Sp+Sh); % herbivore
Mi=y(Sp+Sh+1:Sp+Sh+Sm); % pollinator
%%
% Interaction relationship matrix
matrix_PH=logical(alpha_PH);
matrix_PM=logical(gamma_PM);
%%
% The antagonistic effect of herbivores on plants
F_ph=matrix_PH.*(Pi*Hi')./(1+h*repmat((matrix_PH'*Pi)',[Sp,1]));

% The mutualistic effect of pollinator that acquisition from the plants
E1_pm=matrix_PM.*(Pi*Mi')./(1+h*repmat((matrix_PM'*Pi)',[Sp,1]));

% The mutualistic effects of plant that acquisition from the pollinator
E2_mp=matrix_PM'.*(Mi*Pi')./(1+h*repmat((matrix_PM*Mi)',[Sm,1]));
%%
% Intrinsic growth rate of plants, herbivore and pollinator
r_P=r(1:Sp); 
r_H=r(Sp+1:Sp+Sh);
r_M=r(Sp+Sh+1:Sp+Sh+Sm);
%%
% Expressions of differential equations
dPi=Pi.*(r_P-beta_P*Pi)+(sum(gamma_PM'.*E2_mp))'-sum(alpha_PH.*F_ph,2);
dHi=Hi.*(r_H-beta_H*Hi)+(sum(e*alpha_PH.*F_ph))';
dMi=Mi.*(r_M-beta_M*Mi)+(sum(gamma_PM.*E1_pm))';

dydt(1:Sp)=dPi;
dydt(Sp+1:Sp+Sh)=dHi;
dydt(Sp+Sh+1:Sp+Sh+Sm)=dMi;

% Temporarily protect species from reducing their biomass to too low.
% The trigger probability is very small.
for i=1:Sp+Sh+Sm
    if y(i)<1e-15 && dydt(i)<0
        dydt(i)=0;
    end
end