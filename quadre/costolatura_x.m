function out=costolatura_x(ndlist,planes,y1,y2,z1,z2)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
for ik=1:numel(planes)
    n_c_y1=ndlist(ndlist(find((ndlist(:,2)==planes(ik)).*(ndlist(:,3)==y1)==1)),[1 4]);    
    n_c_z1=ndlist(ndlist(find((ndlist(:,2)==planes(ik)).*(ndlist(:,4)==z1).*(ndlist(:,3)<=y2)==1)),[1 3]);
    n_c_z2=ndlist(ndlist(find((ndlist(:,2)==planes(ik)).*(ndlist(:,4)==z2).*(ndlist(:,3)<=y2)==1)),[1 3]);

    [sortay,i_sort_y1]=sort(n_c_y1(:,2),'ascend');
    [sortaz1,i_sort_z1]=sort(n_c_z1(:,2),'descend');
    [sortaz2,i_sort_z2]=sort(n_c_z2(:,2),'ascend');
    
    ny1=numel(n_c_y1(:,1));
    nz1=numel(n_c_z1(:,1));
    nz2=numel(n_c_z2(:,1));
           
    node_beam_cerc=[n_c_z1(i_sort_z1,1)' n_c_y1(i_sort_y1,1)' n_c_z2(i_sort_z2,1)'];
    beam_orient=[180*ones(1,nz1) 90*ones(1,ny1) 0*ones(1,nz2)];
    
    
    iz=1;
    for ij=1:numel(node_beam_cerc)-1

        if node_beam_cerc(ij)~=node_beam_cerc(ij+1)
            beam(iz,:)= [node_beam_cerc(ij) node_beam_cerc(ij+1) beam_orient(ij)];
            iz=iz+1;
        end
    end   
    
    out{ik}=beam;
               
end
