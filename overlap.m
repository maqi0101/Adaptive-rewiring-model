%This function is used to calculate the niche overlap between two species

function v_ol = overlap(s1,s2,nw,trSpan)

% To avoid edge effects, the niche axis is defined circular (''periodic''),
% so that each species has equal numbers of partners or competitors on both sides
% Refer to Weiran Cai(2020); Marten Scheffer(2006)

    v_ol = exp(-min([abs(s1-s2),trSpan-abs(s1-s2)])^2/(4*nw^2));
    
end