function out=costolatura_z(ndlist,planes,y1,y2,x1,x2)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
for ik=1:numel(planes)
    n_c_y1=ndlist(ndlist(find((ndlist(:,4)==planes(ik)).*(ndlist(:,3)==y1)==1)),[1 2]);    
    n_c_x1=ndlist(ndlist(find((ndlist(:,4)==planes(ik)).*(ndlist(:,2)==x1).*(ndlist(:,3)<=y2)==1)),[1 3]);
    n_c_x2=ndlist(ndlist(find((ndlist(:,4)==planes(ik)).*(ndlist(:,2)==x2).*(ndlist(:,3)<=y2)==1)),[1 3]);

    [sortay1,i_sort_y1]=sort(n_c_y1(:,2),'ascend');
    [sortax1,i_sort_x1]=sort(n_c_x1(:,2),'descend');
    [sortax2,i_sort_x2]=sort(n_c_x2(:,2),'ascend');
    
    ny1=numel(n_c_y1(:,1));
    nx1=numel(n_c_x1(:,1));
    nx2=numel(n_c_x2(:,1));
    
    
    
    node_beam_cerc=[n_c_x1(i_sort_x1,1)' n_c_y1(i_sort_y1,1)' n_c_x2(i_sort_x2,1)'];
    beam_orient=[90*ones(1,nx1) 90*ones(1,ny1) 90*ones(1,nx2)];
    
    
    iz=1;
    for ij=1:numel(node_beam_cerc)-1

        if node_beam_cerc(ij)~=node_beam_cerc(ij+1)
            beam(iz,:)= [node_beam_cerc(ij) node_beam_cerc(ij+1) beam_orient(ij)];
            iz=iz+1;
        end
    end   
    
    out{ik}=beam;
               
end
%%
%%
