function newndlist=roundtoalpharef(ndlist,alpharef)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
for ij=1:numel(ndlist(:,1))
    alp=atan2(ndlist(ij,3),ndlist(ij,2));
    R_temp=sqrt(ndlist(ij,3)^2+ndlist(ij,2)^2);    
        if alp~=0      
            
            alp=alpharef*round(alp/(alpharef));
        end
    newndlist(ij,:)=[ndlist(ij,1),R_temp*cos(alp),R_temp*sin(alp),ndlist(ij,4)];
end
