function polarquad=polarnode(circular_mesh)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
for ij=1:numel(circular_mesh(:,1))
    R=round(sqrt(circular_mesh(ij,2)^2+circular_mesh(ij,3)^2));
    alp=atan2(circular_mesh(ij,3),circular_mesh(ij,2));
    if alp<0
        alp=2*pi+alp;
    end
    polarlist(ij,:)=[circular_mesh(ij,1) R alp];
end


   
polarlist(:,3)=round(polarlist(:,3)*1000)/1000;
alplist=unique(polarlist(:,3));
alplist=sort(alplist);
polarquad=[];
%numel(alplist)-1;
for ij=1:numel(alplist)-1;
    polarselecta_1=polarlist(polarlist(:,3)==alplist(ij),:);
    polarselecta_2=polarlist(polarlist(:,3)==alplist(ij+1),:);
    polarselecta=[polarselecta_1;polarselecta_2];
    polarquad=[polarquad; quadronode2(polarselecta(:,[1 3 2]))];
end

polarselecta_1=polarlist(polarlist(:,3)==alplist(1),:);
polarselecta_2=polarlist(polarlist(:,3)==alplist(end),:);
polarselecta=[polarselecta_1;polarselecta_2];
lastsector=quadronode2(polarselecta(:,[1 3 2]));
lastsector=lastsector(:,[4 3 2 1 5]);
polarquad=[polarquad;lastsector ];
end
