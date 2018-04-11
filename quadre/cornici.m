function out=cornici(cornice,maxX,minX,maxZ,minZ)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------

cornix0=cornice(cornice(:,2)<=minX,[1 2 4]);
cor_x0=quadronode(cornix0);
cor_x0=[ones(numel(cor_x0(:,1)),1) 3*ones(numel(cor_x0(:,1)),1) cor_x0];

cornixL=cornice(cornice(:,2)>=maxX,[1 2 4]);
cor_xL=quadronode(cornixL);
cor_xL=[ones(numel(cor_xL(:,1)),1) 3*ones(numel(cor_xL(:,1)),1) cor_xL];

corniz0=cornice(cornice(:,4)<=minZ,[1 2 4]);
cor_z0=quadronode(corniz0);
cor_z0=[ones(numel(cor_z0(:,1)),1) 3*ones(numel(cor_z0(:,1)),1) cor_z0];
% 
cornizL=cornice(cornice(:,4)>=maxZ,[1 2 4]);
cor_zL=quadronode(cornizL);
cor_zL=[ones(numel(cor_zL(:,1)),1) 3*ones(numel(cor_zL(:,1)),1) cor_zL];

out=[cor_x0; cor_xL; cor_z0; cor_zL]; 
out=unique(out,'rows');
