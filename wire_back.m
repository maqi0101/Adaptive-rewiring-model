% This function is used to restore the interaction

function [alpha_PH,gamma_PM] = wire_back(alpha_PH,gamma_PM,int_old,idstart,idend_new,idend_old,Type)

if Type==1 % if this rewiring is antagonistic
    alpha_PH(idend_new,idstart)=0; % disconnect the new interaction  
    alpha_PH(idend_old,idstart)=int_old; % restore the original interaction
    
else % if this rewiring is mutualistic
    gamma_PM(idend_new,idstart)=0; 
    gamma_PM(idend_old,idstart)=int_old; 
end 