function quad=quadronode(selecta)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
[sorta,isorta]=sort(selecta(:,2));
selecta=selecta(isorta,:);
n_el_1=numel(find(selecta(:,2)==selecta(1,2)))-1;
% n_el_2=numel(find(selecta(:,3)==selecta(1,3)))-1;
n_el_2=numel(selecta(:,1))/(n_el_1+1)-1;

i_nod=reshape(selecta(:,1),n_el_1+1,n_el_2+1);
% dir1=reshape(selecta(:,2),n_el_1+1,n_el_2+1);
dir2=reshape(selecta(:,3),n_el_1+1,n_el_2+1);

for ij=1:numel(dir2(1,:))
    [sorta2,sortdir2]=sort(dir2(:,ij));
    i_nod(:,ij)=i_nod(sortdir2,ij);
    dir2(:,ij)=dir2(sortdir2,ij);
end
 
iz=1;
for ij=1:numel(i_nod(:,1))-1
    for ik=1:numel(i_nod(1,:))-1
        % nodo1 nodo2 nodo3 nodo4 ordinata media
        quad(iz,:)=[i_nod(ij,ik) i_nod(ij+1,ik) i_nod(ij+1,ik+1) i_nod(ij,ik+1) (dir2(ij,ik)+dir2(ij+1,ik+1))/2];
        iz=iz+1;
    end
end
end