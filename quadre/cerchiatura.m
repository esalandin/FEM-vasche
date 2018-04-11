function out=cerchiatura(ndlist,planes,x1,x2,z1,z2)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
for ik=1:numel(planes)
    n_c_x1=ndlist(find((ndlist(:,3)==planes(ik)).*(ndlist(:,2)==x1)==1),[1 4]);
    n_c_x2=ndlist(find((ndlist(:,3)==planes(ik)).*(ndlist(:,2)==x2)==1),[1 4]);
    n_c_z1=ndlist(find((ndlist(:,3)==planes(ik)).*(ndlist(:,4)==z1)==1),[1 2]);
    n_c_z2=ndlist(find((ndlist(:,3)==planes(ik)).*(ndlist(:,4)==z2)==1),[1 2]);
    nx1=numel(n_c_x1(:,1));
    nx2=numel(n_c_x2(:,1));
    nz1=numel(n_c_z1(:,1));
    nz2=numel(n_c_z2(:,1));

    [not_used_1,i_sort_x1]=sort(n_c_x1(:,2),'ascend');
    [not_used_2,i_sort_x2]=sort(n_c_x2(:,2),'descend');
    [not_used_3,i_sort_z1]=sort(n_c_z1(:,2),'descend');
    [not_used_4,i_sort_z2]=sort(n_c_z2(:,2),'ascend');

    node_beam_cerc=[n_c_x1(i_sort_x1,1)' n_c_z2(i_sort_z2,1)' n_c_x2(i_sort_x2,1)' n_c_z1(i_sort_z1,1)' ];
    beam_orient=[0*ones(1,nx1) 0*ones(1,nz2) 0*ones(1,nx2) 180*ones(1,nz1-1)];

    
    
    iz=1;
    for ij=1:numel(node_beam_cerc)-1

        if node_beam_cerc(ij)~=node_beam_cerc(ij+1)
            beam(iz,:)= [node_beam_cerc(ij) node_beam_cerc(ij+1) beam_orient(ij)];
            iz=iz+1;
        end
    end   
    
    out{ik}=beam;
               
end
