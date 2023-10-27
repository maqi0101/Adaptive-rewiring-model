% This function is used to rewiring one link

function [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old,int_old]=rewiring(Type,N,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta)

% rewiring one link

if (Type==1)  % rewiring H->P: herbivore to plant
    idH=randi(Sh); % select a species randomly from herbivores
    while isempty(find(alpha_PH(:,idH),1))  % if idH has no link to plant at all£¬reselect
        idH=randi(Sh);
    end
    
    gain_old=N(Sp+idH);  % the biomass of this animal before rewiring
    
    idP_nz=find(alpha_PH(:,idH));  % find out all of the animal's plant partners
    idP_old=idP_nz(randi(length(idP_nz))); % randomly select a plant from it
    
    int_old=alpha_PH(idP_old,idH); % the strength of the interaction between them
    
    idP_new=idP_old;
    deg_old=nnz(alpha_PH(idP_old,:)); % the number of antagonistic partners of this plant
    pr=1/deg_old^eta;  % rewiring probability
    if rand>pr
        idP_z=find(alpha_PH(:,idH)==0); % find out all plants which  the animal unconnected
        if ~isempty(idP_z)
            idP_new=idP_z(randi(length(idP_z))); % randomly select a plant from it
            
            % rewire to a new plant
            alpha_PH(idP_old,idH)=0; % disconnect the original interaction
            alpha_PH(idP_new,idH)=int_ant*overlap(trait(idP_new),trait(Sp+idH),nw,trSpan); % rewiring
        end
    end
    
    idstart=idH; % record the number of the animal
    idend_old=idP_old; % record the number of disconnected plants
    idend_new=idP_new; % record the number of new plants
    
else  % rewiring M->P: pollinator to  plant
    % The following steps are the same as before
    idM=randi(Sm);
    while isempty(find(gamma_PM(:,idM),1))
        idM=randi(Sm);
    end
    
    gain_old=N(Sp+Sh+idM);
    
    idP_nz=find(gamma_PM(:,idM));
    idP_old=idP_nz(randi(length(idP_nz)));
    
    int_old=gamma_PM(idP_old,idM);
    
    idP_new=idP_old;
    deg_old=nnz(gamma_PM(idP_old,:));
    pr=1/deg_old^eta;
    
    if rand>pr
        idP_z=find(gamma_PM(:,idM)==0);
        if ~isempty(idP_z)
            idP_new=idP_z(randi(length(idP_z)));
            
            % rewire to a new plant
            gamma_PM(idP_old,idM)=0;
            gamma_PM(idP_new,idM)=int_mut*overlap(trait(idP_new),trait(Sp+Sh+idM),nw,trSpan);
        end
    end
    
    idstart=idM;
    idend_old=idP_old;
    idend_new=idP_new;
    
end

