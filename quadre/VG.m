function VG(nome,L_x,L_y,L_z,s_f,s_p,sb,B,D,tt,c_x,c_y,c_z,mat,path,descr,rho_l,apg_type,apg_pos,ad_nd_x,ad_nd_z,dx,dy,dz)

% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
% File per realizzazione del modello di una vasca galvanica date le
% dimensioni L (lunghezza), H (altezza), W (profondità), s (spessore del 
% fondo e delle pareti della vasca), sb (spessore del bordo), Lb (lunghezza
% del bordo)
% global keepvalue keepvalue2

% Per compensare il fatto che gli elementi shell vengono creati con i nodi sul piano medio, si "aumentano" 
% fittiziamente le dimensioni della vasca per avere il reale interno vasca
L_x=L_x+s_p;
L_z=L_z+s_p;
L_y=L_y+s_p;

% offset livello di liquido da bordo vasca
offset_h2o=0;
y_h2o=L_y-offset_h2o;

% discretizzazione elementi
%dxdef=0.04;
%dydef=0.04;
%dzdef=0.04;

% proprietà tubolare di rinforzo [m]
T1=tt;         
T2=tt;


% parametri bordo
% si ipotizza che il bordo occupi lo stesso spazio del tubolare di rinforzo
% della cerchiatura
Lb=B;
Hb=0.03;

% Rb=0.2;

%% Scelta del materiale che compone la vasca

% Scelta da parte dell'utente del materiale da impiegare per la
% realizzazione della vasca galvanica, la scelta è fra PP, acciaio e
% Titanio. Tutta la vasca viene realizzata con lo stesso materiale
% mat=input('Scegli il materiale che costituisce la vasca (PP=1 Fe360=2 Ti=3):  ');

% mat=2;

if mat==1
    E=1600;      %[MPa]
    rho=900;     %[kg/m^3]
    ni=0.3;
    spec_heat=519.163;   %[J/kg/K]
    cond=0.22;   %[W/m*K]  
    materiale='PP';
% Lo spessore del bordo per le vasche in PP viene sovrascritto
    sb=0.04; 
elseif mat==2
    E=198000;   %[MPa]
    rho=7870;    %[kg/m^3]
    ni=0.3;
    spec_heat=500;   %[J/kg/K]
    cond=16;   %[W/m*K] 
    materiale='Acciaio';
elseif mat==3
    E=103000;   %[MPa]
    rho=4510;    %[kg/m^3]
    ni=0.3;
    spec_heat=519.163;   %[J/kg/K]
    cond=21.788;   %[W/m*K]
    materiale='Titanio';
elseif mat==4
    E=1300;   %[MPa]
    rho=950;    %[kg/m^3]
    ni=0.42;
    spec_heat=1900;   %[J/kg/K]
    cond=0.4;   %[W/m*K]
    materiale='PE';    
end


if isnan(c_x)==1
    pos_x=sort([0 L_x]);
else
    pos_x=sort([0 c_x L_x]);
end

if isnan(ad_nd_x)~=1   
    pos_x=sort([pos_x ad_nd_x]);
end

if strcmp(apg_type,'HEA200')==1
    E_tr=198000;   %[MPa]
    rho_tr=7870;    %[kg/m^3]
    ni_tr=0.3;
    spec_heat_tr=500;   %[J/kg/K]
    cond_tr=16;   %[W/m*K] 
%     materiale_tr='Acciaio';        
    pos_x=sort([pos_x apg_pos apg_pos+0.1 apg_pos-0.1]);    
end

if strcmp(apg_type,'HEA100')==1
    E_tr=198000;   %[MPa]
    rho_tr=7870;    %[kg/m^3]
    ni_tr=0.3;
    spec_heat_tr=500;   %[J/kg/K]
    cond_tr=16;   %[W/m*K] 
%     materiale_tr='Acciaio';        
    pos_x=sort([pos_x apg_pos apg_pos+0.05 apg_pos-0.05]);    
end

x=nodes(pos_x,dx);

if isnan(c_y)==1
    pos_y=sort([0 L_y]);
else
    pos_y=sort([0 c_y L_y]);
end

y=nodes(pos_y,dy);

if isnan(c_z)==1
    pos_z=sort([0 L_z]);
else
    pos_z=sort([0 c_z L_z]);
end

if isnan(ad_nd_z)~=1   
    pos_z=sort([pos_z ad_nd_z]);
end

z=nodes(pos_z,dz);



[nodes_x , nodes_y]=meshgrid(x,y);
nodes_xy=[nodes_x(:) nodes_y(:)];
% break
nodes_xy=nodes_xy((nodes_xy(:,1)~=0).*(nodes_xy(:,1)~=L_x)==1,:);   %elimina nodi di spigolo laterale xy

[nodes_z , nodes_y]=meshgrid(z,y);
nodes_zy=[nodes_z(:) nodes_y(:)];

[nodes_x , nodes_z]=meshgrid(x,z);
nodes_xz=[nodes_x(:) nodes_z(:)];
logical_xz=(nodes_xz(:,1)~=0).*(nodes_xz(:,1)~=L_x).*(nodes_xz(:,2)~=0).*(nodes_xz(:,2)~=L_z);   %elimina tutti nodi spigolo base
nodes_xz=nodes_xz(logical_xz==1,:);



% fondo
ndlist1=[nodes_xz(:,1) zeros(numel(nodes_xz(:,1)),1) nodes_xz(:,2)];   %inserisce colonna zeri tra coord x e z
%pareti xy
ndlist2=[nodes_xy(:,1) nodes_xy(:,2) zeros(numel(nodes_xy(:,1)),1)];    %inserisce colonna zeri dopo coord x e y
ndlist3=[nodes_xy(:,1) nodes_xy(:,2) L_z*ones(numel(nodes_xy(:,1)),1)]; %inserisce colonna di Lz dopo coord x e y
%pareti zy
ndlist4=[zeros(numel(nodes_zy(:,1)),1) nodes_zy(:,2) nodes_zy(:,1)];     %inserisce colonna di zeri prima coord z e y
ndlist5=[L_x*ones(numel(nodes_zy(:,1)),1) nodes_zy(:,2) nodes_zy(:,1)];  %inserisce colonna di Lx prima coord z e y

ndlist=[ndlist1;ndlist2;ndlist3;ndlist4;ndlist5];
ndlist=[(1:numel(ndlist(:,1)))' ndlist];      %aggiunge ID NODI


% elementi fondo vasca
selecta=ndlist(ndlist(:,3)==0,:);
selecta=selecta(:,[1 2 4]);
quad_fondo=quadronode(selecta);
% gruppo proprietà nodi(4) coordinata y media
quad_fondo=[ones(numel(quad_fondo(:,1)),1) ones(numel(quad_fondo(:,1)),1) quad_fondo];
quad_fondo(:,7)=zeros(numel(quad_fondo(:,7)),1);

% elementi parete x=0
selecta=ndlist(ndlist(:,2)==0,:);
selecta=selecta(:,[1 4 3]);
quad_x0=quadronode(selecta);
quad_x0=[ones(numel(quad_x0(:,1)),1) 2*ones(numel(quad_x0(:,1)),1) quad_x0];

% elementi parete x=L_x
selecta=ndlist(ndlist(:,2)==L_x,:);
selecta=selecta(:,[1 4 3]);
% nodo1 nodo2 nodo3 nodo4 ordinata media
quad_lx=quadronode(selecta);
% gruppo proprietà nodo1 nodo2 nodo3 nodo4 ordinata media
quad_lx=[ones(numel(quad_lx(:,1)),1) 2*ones(numel(quad_lx(:,1)),1) quad_lx(:,4:-1:1) quad_lx(:,5)];

% elementi parete z=0
selecta=ndlist(ndlist(:,4)==0,:);
selecta=selecta(:,[1 2 3]);
quad_z0=quadronode(selecta);
quad_z0=[ones(numel(quad_z0(:,1)),1) 2*ones(numel(quad_z0(:,1)),1) quad_z0(:,4:-1:1) quad_z0(:,5)];

% elementi parete z=L_z
selecta=ndlist(ndlist(:,4)==L_z,:);
selecta=selecta(:,[1 2 3]);
quad_lz=quadronode(selecta);
quad_lz=[ones(numel(quad_lz(:,1)),1) 2*ones(numel(quad_lz(:,1)),1) quad_lz];

% qui dentro ho una riga per ogni elemento quad, con coordinate nodi (4)
% valore ordinata media (1) e all'inizio 1 1 o 1 2
quad=[quad_fondo;quad_x0;quad_lx;quad_z0;quad_lz];   



% bordo
% selecta=ndlist(ndlist(:,3)==L_y,:);
% selecta_top_vasca=selecta;
b_x0=ndlist(find((ndlist(:,3)==L_y).*(ndlist(:,2)==0)==1),2:end);   %trova nodi su bordo superiore x=0
b_x1=ndlist(find((ndlist(:,3)==L_y).*(ndlist(:,2)==L_x)==1),2:end); %trova nodi su bordo superiore x=Lx
b_z0=ndlist(find((ndlist(:,3)==L_y).*(ndlist(:,4)==0)==1),2:end); %trova nodi su bordo superiore z=0
% .*(ndlist(:,2)~=L_x).*(ndlist(:,2)~=0) 
b_z1=ndlist(find((ndlist(:,3)==L_y).*(ndlist(:,4)==L_z)==1),2:end);    %trova nodi su bordo superiore z=Lz
% .*(ndlist(:,2)~=L_x).*(ndlist(:,2)~=0)

D_l=(Lb/(round(Lb/dx))):(Lb/(round(Lb/dx))):Lb;   %calcola lato elementi bordo

b_x0=sort([b_x0;b_x0(1:numel(D_l),1:2) -D_l';b_x0(1:numel(D_l),1:2) max(b_x0(:,3))+D_l']);   %aggiunge nodi bordo a x=0
b_x1=sort([b_x1;b_x1(1:numel(D_l),1:2) -D_l';b_x1(1:numel(D_l),1:2) max(b_x1(:,3))+D_l']);  %aggiunge nodi bordo a x=Lx


% bor_nod=[b_x0; b_x1; b_z0; b_z1];
bor_nod=[];
for ij=1:numel(D_l)
    Dl=D_l(ij);
    bor_nod=[bor_nod;b_x0+repmat(Dl*[-1 0 0],numel(b_x0(:,1)),1)];   %aggiunge nodi fuori bordo a x=0
    bor_nod=[bor_nod;b_x1+repmat(Dl*[1 0 0],numel(b_x1(:,1)),1)];    %aggiunge nodi fuori bordo a x=xL
    bor_nod=[bor_nod;b_z0+repmat(Dl*[0 0 -1],numel(b_z0(:,1)),1)];    %aggiunge nodi fuori bordo a z=0
    bor_nod=[bor_nod;b_z1+repmat(Dl*[0 0 1],numel(b_z1(:,1)),1)];     %aggiunge nodi fuori bordo a z=zL
    PP_bor_nod=bor_nod;   %se vasca di PP no piega bordo
end
bor_nod=[bor_nod; bor_nod-repmat(Hb*[0 1 0],numel(bor_nod(:,1)),1)];  %aggiunge nodi piega su tutto bordo vasca
bor_nod=[[numel(ndlist(:,1))+1:numel(ndlist(:,1))+numel(bor_nod(:,1))]' bor_nod]; %aggiunge ID partendo dagli ID di ndlist
PP_bor_nod=[[numel(ndlist(:,1))+1:numel(ndlist(:,1))+numel(PP_bor_nod(:,1))]' PP_bor_nod]; %aggiunge ID partendo dagli ID di ndlist

% cerchiature NB: i nodi e gli elementi sono gia stati generati per avere
% nodi a quell'y
if isnan(c_y)~=1
    mat_cer=cerchiatura(ndlist,c_y,0,L_x,0,L_z);    
else
    mat_cer=[];
end
% costolature x (simile cerchiature)

if isnan(c_x)~=1
    mat_cosx=costolatura_x(ndlist,c_x,0,max(c_y),0,L_z);
else
    mat_cosx=[];
end
% costolature z
if isnan(c_z)~=1
    mat_cosz=costolatura_z(ndlist,c_z,0,max(c_y),0,L_x);
else
    mat_cosz=[];
end


allbeam=[mat_cer mat_cosx mat_cosz ];   %tutti i nodi travi costolature/cerchiature con relativi angoli per orientare assi principali


if mat==2 || mat==3 % aggiunge nodi cornice e se mat richiede anche cornice verticale
    ndlist=[ndlist; bor_nod];
    cornice_xy1=ndlist(ndlist(:,3)==L_y,:); 
    cor_1=cornici(cornice_xy1,L_x,0,L_z,0); % ID vertici elementi cornice orizzontale

    cor2=[];
    cornice_y=ndlist(ndlist(:,2)==max(ndlist(:,2)),[1 3 4]);
    cor2=[cor2 ; quadronode(cornice_y)]; % ID vertici elementi cornice vert

    cornice_y=ndlist(ndlist(:,2)==min(ndlist(:,2)),[1 3 4]);
    cor2=[cor2 ; quadronode(cornice_y)]; % ID vertici elementi cornice vert

    cornice_y=ndlist(ndlist(:,4)==max(ndlist(:,4)),[1 3 2]);
    cor2=[cor2 ; quadronode(cornice_y)]; % ID vertici elementi cornice vert

    cornice_y=ndlist(ndlist(:,4)==min(ndlist(:,4)),[1 3 2]);
    cor2=[cor2 ; quadronode(cornice_y)]; % ID vertici elementi cornice vert
    cor2=[ones(numel(cor2(:,1)),1) 3*ones(numel(cor2(:,1)),1) cor2];
%             ndlist(ndlist(:,2)==min(ndlist(:,2)),:); 
%             ndlist(ndlist(:,4)==max(ndlist(:,4)),:);
%             ndlist(ndlist(:,4)==min(ndlist(:,4)),:) 
    quad=[quad; cor_1; cor2];   %aggiunge elementi cornice verticale
else   %caso senza cornice verticale
    ndlist=[ndlist; PP_bor_nod];
    cornice_xy1=ndlist(ndlist(:,3)==L_y,:);
    cor_1=cornici(cornice_xy1,L_x,0,L_z,0);    

%     cor2=[];
%     cornice_y=ndlist(ndlist(:,2)==max(ndlist(:,2)),[1 3 4]);
%     cor2=[cor2 ; quadronode(cornice_y)];
% 
%     cornice_y=ndlist(ndlist(:,2)==min(ndlist(:,2)),[1 3 4]);
%     cor2=[cor2 ; quadronode(cornice_y)];
% 
%     cornice_y=ndlist(ndlist(:,4)==max(ndlist(:,4)),[1 3 2]);
%     cor2=[cor2 ; quadronode(cornice_y)];
% 
%     cornice_y=ndlist(ndlist(:,4)==min(ndlist(:,4)),[1 3 2]);
%     cor2=[cor2 ; quadronode(cornice_y)];
%     cor2=[ones(numel(cor2(:,1)),1) 3*ones(numel(cor2(:,1)),1) cor2];
 
    quad=[quad; cor_1];
end

% id_start_bordi=numel(quad(:,1));



start=1;
beam_prop_type=1;

for ij=1:numel(allbeam)
    temp_mat=allbeam{ij};       %selezione matrice costolex, z o cerchiatura
    count=start:start+numel(temp_mat(:,1))-1;   %
    n=numel(count);
%     if mat==1
%         if ij==numel(allbeam)
%             beam_prop_type=2;
%         end
%     end
    temp_mat=[count' ij*ones(n,1) beam_prop_type*ones(n,1) temp_mat];

    allbeam{ij}=temp_mat;
    start=start+numel(temp_mat(:,1));
end


%% inserimento tipo appoggio vasca 
if strcmp(apg_type,'HEA200')==1 || strcmp(apg_type,'HEA100')==1
%  selezione nodi fondo
nodi_f=ndlist(ndlist(:,3)==0,:);    

for ij=1:numel(apg_pos)

%     (ndlist(:,2)>=apg_pos(ij)-0.1)
%     nodi_f_apg=nodi_f(find((nodi_f(:,2)<=apg_pos(ij)+0.1).*(nodi_f(:,2)>=apg_pos(ij)-0.1)==1),:);
    index_nodi_apg=(nodi_f(:,2)<=apg_pos(ij)+0.1).*(nodi_f(:,2)>=apg_pos(ij)-0.1)==1;   %seleziona nodi su trave
    nodi_f_apg=nodi_f(index_nodi_apg,:);
    row_nodi=nodi_f_apg(nodi_f_apg(:,4)==0,:);               %seleziona nodi z=0
    row_nodi_zm=nodi_f_apg(nodi_f_apg(:,4)==L_z,:);      %seleziona nodi a fine z=Lz         
    
    nodi_f_apg=nodi_f_apg(:,[1 2 4]);
    quadtra=quadronode(nodi_f_apg);
    quadtra=[ones(numel(quadtra(:,1)),1) 4*ones(numel(quadtra(:,1)),1) quadtra(:,4:-1:1) quadtra(:,5)];
    quad=[quad; quadtra];    %aggiunge nodi appoggio da vincolare
    
    newz=-dz:-dz:-0.500;        
    [new_nodix , new_nodiz]=meshgrid(row_nodi(:,2),newz);
    newcount=numel(ndlist(:,1))+1:1:(numel(ndlist(:,1))+numel(new_nodix(:)));        
    add_nodes=[newcount' new_nodix(:) zeros(numel(new_nodix(:)),1) new_nodiz(:)];   %crea nodi appoggio fuori vasca
    
    ndlist=[ndlist; add_nodes];         %
    add_nodes=[add_nodes;row_nodi];    
    add_nodes=add_nodes(:,[1 2 4]);
    quadtra_zm=quadronode(add_nodes);
    quadtra_zm=[ones(numel(quadtra_zm(:,1)),1) 4*ones(numel(quadtra_zm(:,1)),1) quadtra_zm(:,4:-1:1) quadtra_zm(:,5)]; 
    %ha generato nodi per appoggio  fuori (z=0) trave e li somma poi agli
    %altri nodi (NB identificativo 1 4)
    quad=[quad; quadtra_zm];
    
    newz=L_z+dz:dz:L_z+0.500;        
    [new_nodix , new_nodiz]=meshgrid(row_nodi_zm(:,2),newz);
    newcount=numel(ndlist(:,1))+1:1:(numel(ndlist(:,1))+numel(new_nodix(:)));        
    add_nodes=[newcount' new_nodix(:) zeros(numel(new_nodix(:)),1) new_nodiz(:)];
    
    ndlist=[ndlist; add_nodes];
    add_nodes=[add_nodes;row_nodi_zm];    
    add_nodes=add_nodes(:,[1 2 4]);
    quadtra_zm=quadronode(add_nodes);
    quadtra_zm=[ones(numel(quadtra_zm(:,1)),1) 4*ones(numel(quadtra_zm(:,1)),1) quadtra_zm(:,4:-1:1) quadtra_zm(:,5)];
    quad=[quad; quadtra_zm];
        %ha generato nodi per appoggio  fuori trave (z=Lz) e li somma poi agli
    %altri nodi (NB identificativo 1 4)


end
end
quad=[[1:numel(quad(:,1))]' quad];

%% Densità del liquido contenuto nella vasca

% Scelta da parte dell'utente della densità del liquido contenuto nella
% vasca, dalla quale dipende la distribuzione di pressioni lungo le pareti
% e il fondo della stessa
% density=input('Densità del liquido contenuto nella vasca [kg/m^3]:  ');
density=rho_l;

%% WARNING

% Warning per far capire all'utente che, nel caso scelga di ralizzare una
% vasca in acciaio ha inserito un valore per la lunghezza del bordo
% inferiore a quello impostato per la flangia (che misura 30x50 [mm])
% if mat==2 && Lb<=0.03
%     disp('WARNING: LA FLANGIA COMPENETRA LA VASCA!!!');
% end

%% Scrittura del file txt

% Apertura del file di testo che verrà poi importato in Straus
if isnan(path)
    nome_modello=strcat(nome,'.txt');
else
   nome_modello=strcat(path,nome,'.txt');
end

pn=fopen(nome_modello,'w');
 
%% MODEL INFORMATION

comment=descr;
% comment=strcat(comment,'dimensioni vasca:',num2str(int32(L_x*1000)),'x',num2str(int32(L_y*1000)),'x',num2str(int32(L_z*1000)));
% comment=strcat(comment,'spessore fondo:',num2str(s_f),'Spessore pareti',num2str(s_p));
% comment=strcat(comment,'elementi di rinforzo',num2str(B*1000),'x',num2str(D*1000),'x',num2str(tt*1000));
% comment=strcat(comment,'materiale vasca:',materiale);

    
% Scrittura dell prime informazioni generali sul file di testo, usate da
% Straus per realizzare poi il modello
fprintf(pn,'%s %s\n','FileFormat','Straus7.2.3.3');
fprintf(pn,'%s %s\n','ModelName','"Modello_vasca_galvanica"');
fprintf(pn,'%s %s\n','Title','""');
fprintf(pn,'%s %s\n','Project','""');
fprintf(pn,'%s %s\n','Author','""');
fprintf(pn,'%s %s\n','Reference','""');
fprintf(pn,'%s %s %s %s %s %s %s \n','Comments','dimensioni vasca:',num2str(int32(L_x*1000)),'x',num2str(int32(L_y*1000)),'x',num2str(int32(L_z*1000)));
fprintf(pn,'%s %s %s \n','spessore fondo:',num2str(s_f),'Spessore pareti',num2str(s_p));
fprintf(pn,'%s %s %s %s %s %s  \n','elementi di rinforzo',num2str(B*1000),'x',num2str(D*1000),'x',num2str(tt*1000));
fprintf(pn,'%s %s  \n','materiale vasca:',materiale);
    
%% UNITS    

% Scrittura delle informazioni per le unità di misura utilizzate da Straus
fprintf(pn,'%s %s\n','LengthUnit','m');
fprintf(pn,'%s %s\n','MassUnit','kg');
fprintf(pn,'%s %s\n','EnergyUnit','J');
fprintf(pn,'%s %s\n','PressureUnit','MPa');
fprintf(pn,'%s %s\n','ForceUnit','N');
fprintf(pn,'%s %s\n','TemperatureUnit','K');

%% GROUPS DEFINITIONS

% Definizione dei gruppi
fprintf(pn,'%s %s %s %s\n','Group','1','16711680','"\\Bordo"');
fprintf(pn,'%s %s %s %s\n','Group','2','3355647','"\\Fondo"');
fprintf(pn,'%s %s %s %s\n','Group','3','3407692','"\\Fianchi"');
fprintf(pn,'%s %s %s %s\n','Group','4','3407846','"\\Model"');

%% FREEDOM CASE DEFINITIONS

% Definizione dei freedom case
fprintf(pn,'%s %s %s %s %s\n','FreedomCase','2','0','1','"\\Freedom Case 1"');
    
%% LOAD CASE DEFINITIONS

% Definizione del load case nel quale viene inserita la gravità e
% successivamente le pressioni normali ad ogni elemento plate
fprintf(pn,'%s %s %s %s\n','LoadCase','1','1','"\\Load Case 1"');
fprintf(pn,'%s %s %s\n', 'Gravity','2','-9.81000000000000E+0');
fprintf(pn,'%s %s\n','LCInclude','3');
    
%% COORDINATE SYSTEM DEFINITIONS

% Definizione del sistema di coordinate principale
fprintf(pn,'%s %s %s %s\n','CoordSys','1','"\\Global XYZ"','GlobalXYZ');

%% NODE COORDINATES

% % Stampa della matrice dei nodi su file txt
for i=1:numel(ndlist(:,1));
% Nell'ordine va posto l'id del nodo e le coordinate dello stesso in X,Y,Z
fprintf(pn,'%s %s %s %s %s\n','Node',num2str(ndlist(i,1)),...
     num2str(ndlist(i,2)),num2str(ndlist(i,3)),num2str(ndlist(i,4)));
end

% % Stampa della matrice degli elementi su file txt
jk=1;
for ij=1:numel(allbeam)
    beam=allbeam{ij};
    for ik=1:numel(beam(:,1))
        
%     id gruppo proprietà nodo 1 nodo 2   
        fprintf(pn,'%s %s %s %s %s %s \n','Beam',num2str(beam(ik,1)),num2str(beam(ik,2)),num2str(beam(ik,3)),...
            num2str(beam(ik,4)),num2str(beam(ik,5)));
        jk=jk+1;
    end
end

% rotazione sdr principali beams
for ij=1:numel(allbeam)
    beam=allbeam{ij};
    for ik=1:numel(beam(:,1))        
        fprintf(pn,'%s %s %s  \n','BmAngle',num2str(beam(ik,1)),num2str(beam(ik,end)));
        jk=jk+1;
    end
end

% offset beams

for ij=1:numel(allbeam)
    beam=allbeam{ij};
%     if mat==1
%     if ij==numel(allbeam)
%        offset1=-(B_b/2-0.005); 
%        offset2=-D_b/2;        
%     else
%         offset1=-B/2;
%         offset2=0;
%     end    
%     else
        offset1=-B/2-s_p/2;
        offset2=0;       
%     end
    for ik=1:numel(beam(:,1))        
        fprintf(pn,'%s %s %s %s %s\n','BmOffset',num2str(beam(ik,1)),num2str(offset1),'  ',num2str(offset2));
        jk=jk+1;
    end
end



for i=1:numel(quad(:,1));
    % Nell'ordine va posto l'id del nodo, il gruppo a cui appartiene  
    % l'elemento, la proprietà e i 4 nodi che lo compongono
    fprintf(pn,'%s %s %s %s %s %s %s %s\n','Quad4',...
            num2str(quad(i,1)),num2str(quad(i,2)),num2str(quad(i,3)),...
            num2str(quad(i,4)),num2str(quad(i,5)),num2str(quad(i,6)),num2str(quad(i,7)));
end

%% PLATE OFFSET

for i=1:numel(quad(:,1));
if quad(i,3)==3 || quad(i,3)==4
    if quad(i,3)==3
    PlOf=-0.00;                            %%%!!
    end
    if quad(i,3)==4
    PlOf=0.196/2+s_f/2;
    end
    fprintf(pn,'%s %s %s \n','PlOffset',...
            num2str(quad(i,1)),num2str(PlOf));
end
end


keepvalue=quad;
for i=1:numel(quad(:,1));
%     % Nell'ordine va posto l'id dell'elemento e il valore di pressione
%     % corrispondente
    if quad(i,3)==3 || quad(i,3)==4
        
    else
    p=-9.806*density*(y_h2o-quad(i,8))*10^-6;
        fprintf(pn,'%s %s %s %s %s\n','Load Case 1','PlPressure',...
                '1',num2str(quad(i,1)),...
                num2str(p));
    end
end
% 

% 


% Stampa su file di testo delle varie proprietà necessarie a definire il
% modello della vasca

%% BEAM PROPERTIES

    A=B*D-(B-2*T1)*(D-2*T1);
    I11=1/12*(B*D^3)-1/12*((B-2*T1)*(D-2*T2)^3);
    I22=1/12*(B^3*D)-1/12*((B-2*T1)^3*(D-2*T2));
    J=2*1/3*B*T1^3+2*1/3*D*T2^3;


    fprintf(pn,'%s %s\n','BeamProp','1','"\\Rinforzi"');
    fprintf(pn,'%s %s\n','MaterialName','"Acciaio"');
    fprintf(pn,'%s %s\n','Modulus',num2str(210000));
    fprintf(pn,'%s %s\n','ShearMod',num2str(76000));
    fprintf(pn,'%s %s\n','Poisson',num2str(0.3));
    fprintf(pn,'%s %s\n','UsePoisson','TRUE');
    fprintf(pn,'%s %s\n','Density',num2str(7900));
    fprintf(pn,'%s %s\n','Expansion',num2str(10));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(10));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(10));
    fprintf(pn,'%s %s\n','Area',num2str(A));    
    fprintf(pn,'%s %s\n','MomentI11',num2str(I11));
    fprintf(pn,'%s %s\n','MomentI22',num2str(I22));
    fprintf(pn,'%s %s\n','MomentJ',num2str(1.691800714286e-7));
    fprintf(pn,'%s %s\n','SectionType','HollowRect');
    fprintf(pn,'%s %s\n','B',num2str(B));
    fprintf(pn,'%s %s\n','D',num2str(D));
    fprintf(pn,'%s %s\n','T1',num2str(T1));
    fprintf(pn,'%s %s\n','T2',num2str(T2));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');


%% PLATE PROPERTIES

% Proprietà per gli elementi del fondo
    fprintf(pn,'%s %s %s\n','PlateShellProp','1','"\\Fondo"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni));
    fprintf(pn,'%s %s\n','Density',num2str(rho));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_f));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_f));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');
    
% Proprietà per gli elementi delle pareti
    fprintf(pn,'%s %s %s\n','PlateShellProp','2','"\\Pareti"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni));
    fprintf(pn,'%s %s\n','Density',num2str(rho));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_p));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_p));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');    
 
% Proprietà per gli elementi del bordo    
    fprintf(pn,'%s %s %s\n','PlateShellProp','3','"\\Bordo"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni));
    fprintf(pn,'%s %s\n','Density',num2str(rho));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat));    
    fprintf(pn,'%s %s\n','MemThick',num2str(sb));
    fprintf(pn,'%s %s\n','BendThick',num2str(sb));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');   
    
 if strcmp(apg_type,'HEA200')==1  || strcmp(apg_type,'HEA100')==1  
     IxxHEA200=1/12*(0.2*0.196^3)-1/12*2*((0.2-0.005)/2*(0.196-2*0.01)^3);
     st_eq=(IxxHEA200*3/0.2)^(1/3);  %la trave di appoggio è idealizzata come un plate di inerzia equivalente
% Proprietà per gli elementi del della trave di appoggio sottovasca    
    fprintf(pn,'%s %s %s\n','PlateShellProp','4','"\\Trave"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_tr));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_tr));
    fprintf(pn,'%s %s\n','Density',num2str(rho_tr));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_tr));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_tr));    
    fprintf(pn,'%s %s\n','MemThick',num2str(st_eq));
    fprintf(pn,'%s %s\n','BendThick',num2str(st_eq));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');       
 end

% % Chiusura del file di testo
fclose all;    
disp('vasca creata');
end

