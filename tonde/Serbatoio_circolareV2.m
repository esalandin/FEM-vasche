%--------------------------------------------------------------------------
%   CREAZIONE SERBATOIO CIRCOLARE v. 2
%   ultima modifica 25/09/2017
%   CHANGELOG:
%   - fissato angolo di 15° per la realizzazione della mesh
%   - suddivisione in settori per carico vento
%--------------------------------------------------------------------------

clear all
close all
clc
fclose all;

%% parametri input

% D0 diametro del Serbatoio [mm]
% H altezza del serbatoio [mm]
D0=2500;
H=3500;
% s_f spessore fondo
s_f=30;
% s_p spessore pareti
s_p=30;
% s_t spessore pareti
s_t=30;

% definizione inspessimento in Polipropilene
D=30;
% 
B=150;


% D dimensione maggiore tubolare rinforzo
Dr=150;
% B dimensione minore tubolare rinforzo
Br=40;
% T1 spessore tubolare di rinforzo
T1r=3;


% densità liquido nel serbatoio [kg/m3]
density=1000;
% livello liquido nel serbatoio [mm]
y_h2o=3250;
% condizioni operative del serbatoio [mbar]
% ATTENZIONE: VALORI DI RIFERIMENTO DA VERIFICARE CASO PER CASO
% pressione aggiuntiva a serbatoio pieno [mbar]
pr_max=25;
% depressione a serbatoio vuoto [mbar]
pr_min=-10;
% pressione anteriore del vento [N/m2]
pwind_front=-464;
% depressione posteriore del vento [N/m2]
pwind_back=-199;

% posizioni degli inspessimenti di serbatoio - vettore [mm]
y_c=[1500 2269 3000];

% posizioni degli cerchiature di serbatoio - vettore [mm]
y_cr=[1000 1500 2000];

% rinforzi radiali serbatoio - vettore [mm]
alp0=45;
alp0=alp0/180*pi;
D_alp=60;
D_alp=D_alp/180*pi;




% nome e percorso file da importare in Straus
nome_modello='C:\Users\test\Desktop\serbatoio.txt';

%% fine parametri input

% conversione in [m]
D0=D0*1e-3;
H=H*1e-3;
y_h2o=y_h2o*1e-3;
s_f=s_f*1e-3;
s_p=s_p*1e-3;
s_t=s_t*1e-3;
D=D*1e-3;
B=B*1e-3;
Dr=Dr*1e-3;
Br=Br*1e-3;
T1r=T1r*1e-3;
T2r=T1r;
y_c=y_c*1e-3;
y_cr=y_cr*1e-3;

% conversione in Mpa
pr_max=pr_max*1e-3*0.1;
pr_min=pr_min*1e-3*0.1;
pwind_back=pwind_back*1e-6;
pwind_front=pwind_front*1e-6;

R0=D0/2;

% Lunghezza di riferimento per elementi della mesh [mm]
% alpha_ref angolo al centro
alpha_ref=15; % °
% conversione in radianti
alpha_ref=alpha_ref*pi/180;

% il rinforzo radiale si adatta sulla mesh esistente in base all'angolo di
% discretizzazione scelto
n_alp0=round(alp0/alpha_ref);
n_D_alp=round(D_alp/alpha_ref);
alp0=n_alp0*alpha_ref;
D_alp=n_D_alp*alpha_ref;
alpR=(alp0:D_alp:2*pi+alp0);
alpR=alpR-pi;
morethan2pi=find(alpR>2*pi);
alpR(morethan2pi)=alpR(morethan2pi)-2*pi;


alpha_ref_list=-pi/4:alpha_ref:pi/4;
add_alpha_1=(alpha_ref_list(1)+alpha_ref_list(2) )/2;
add_alpha_2=(alpha_ref_list(end-1)+alpha_ref_list(end) )/2;
alpha_ref_list=[alpha_ref_list add_alpha_2 add_alpha_1 0];
alpha_ref_list=unique(alpha_ref_list);
alpha_ref_list=sort(alpha_ref_list);
L_el=(D0/2)/((pi/2)/alpha_ref);
% L_el=50;
% L_el=L_el*1e-3;

%Proprietà polipropilene
    E_pp=1600;      %[MPa]
    rho_pp=900;     %[kg/m^3]
    ni_pp=0.3;
    spec_heat_pp=519.163;   %[J/kg/K]
    cond_pp=0.22;   %[W/m*K]  
    materiale_pp='PP';
    
%Proprietà acciaio
    E=198000;   %[MPa]
    rho=7870;    %[kg/m^3]
    ni=0.3;
    spec_heat=500;   %[J/kg/K]
    cond=16;   %[W/m*K] 
    materiale='Acciaio';           


%% mesh centrale (quadrata)

% [x_fondo1,y_fondo1]=meshgrid(-round(D0/(4*L_el))*L_el:L_el:round(D0/(4*L_el))*L_el,-round(D0/(4*L_el))*L_el:L_el:round(D0/(4*L_el))*L_el);

[x_fondo,y_fondo]=meshgrid(D0/4*tan(alpha_ref_list),D0/4*tan(alpha_ref_list));
griglia_fondo=[x_fondo(:) y_fondo(:)];
z_fondo=zeros(numel(griglia_fondo(:,1)),1);
griglia_fondo=[griglia_fondo z_fondo];

%% mesh periferica

% nodi su circonferenza
columns=find(abs(griglia_fondo(:,1))==max(griglia_fondo(:,1))==1);
rows=find(abs(griglia_fondo(:,2))==max(griglia_fondo(:,2))==1);
index_cornice=[columns rows];
index_cornice=unique(index_cornice);
griglia_cornice=griglia_fondo(index_cornice,:);
ic=sort(index_cornice);

ii=1;
jj=1;
for ij=1:numel(griglia_fondo(:,1))
    if ij==ic(ii)
        ii=ii+1;
    else
        griglia_interna(jj,:)=griglia_fondo(ij,:);
        jj=jj+1;
    end
end

Rstart=(max(griglia_fondo(:,1))^2+max(griglia_fondo(:,2))^2)^(1/2)+L_el;
ii=1;
n_elc=round((R0-D0/4)/L_el);
for ij=1:numel(griglia_cornice(:,1))
   alp=atan2(griglia_cornice(ij,2),griglia_cornice(ij,1));

   keepalp(ij)=alp;

   R1=sqrt(griglia_cornice(ij,2)^2+griglia_cornice(ij,1)^2);
   
   L_elc=(R0-R1)/n_elc;
   R_c=R1+L_elc:L_elc:R0-L_elc;
%    numel(R_c)
%    R_c(end)=R0;
   griglia_bordo(ij,:)=[R0*cos(alp) R0*sin(alp) 0];
   
   for ik=1:numel(R_c)        
        griglia_esterna(ii,:)=[R_c(ik)*cos(alp) R_c(ik)*sin(alp) 0];       
       ii=ii+1;
   end
end

% nodi fondo serbatoio
% ndlist1 griglia interna (quadrata) senza la cornice esterna
ndlist1=[[1:numel(griglia_interna(:,1))]' griglia_interna];
% ndlist2 cornice della griglia interna (quadrata)
ndlist2=[[ndlist1(end,1)+1:ndlist1(end,1)+numel(griglia_cornice(:,1))]' griglia_cornice];
% ndlist3 griglia esterna (transizione circolare quadrata)
ndlist3=[[ndlist2(end,1)+1:ndlist2(end,1)+numel(griglia_esterna(:,1))]' griglia_esterna];
% ndlist1 nodi esterni fondo circolare r=R0
ndlist4=[[ndlist3(end,1)+1:ndlist3(end,1)+numel(griglia_bordo(:,1))]' griglia_bordo];
% nodi top serbatoio
% ndlist5 griglia interna (quadrata) senza la cornice esterna TOP
ndlist5=[[ndlist4(end,1)+1:ndlist4(end,1)+numel(griglia_interna(:,1))]' griglia_interna];
% ndlist6 cornice della griglia interna (quadrata) TOP
ndlist6=[[ndlist5(end,1)+1:ndlist5(end,1)+numel(griglia_cornice(:,1))]' griglia_cornice];
% ndlist7 griglia esterna (transizione circolare quadrata) TOP
ndlist7=[[ndlist6(end,1)+1:ndlist6(end,1)+numel(griglia_esterna(:,1))]' griglia_esterna];
% ndlist8 nodi esterni fondo circolare r=R0
ndlist8=[[ndlist7(end,1)+1:ndlist7(end,1)+numel(griglia_bordo(:,1))]' griglia_bordo];
% fisso z=H per tutti i nodi TOP
ndlist5(:,4)=H*ones(numel(ndlist5(:,1)),1);
ndlist6(:,4)=H*ones(numel(ndlist6(:,1)),1);
ndlist7(:,4)=H*ones(numel(ndlist7(:,1)),1);
ndlist8(:,4)=H*ones(numel(ndlist8(:,1)),1);
% nodi pareti serbatoio

coord_z=[L_el H-L_el y_c y_cr];
coord_z=sort(coord_z);
n_z=[];
for ij=1:numel(coord_z)-1
    n_el=round((coord_z(ij+1)-coord_z(ij))/L_el);
    if n_el==0
        n_el=1;
    end
    seq_z=coord_z(ij):(coord_z(ij+1)-coord_z(ij))/n_el:coord_z(ij+1);
    n_z=[n_z seq_z];    
end
n_z=unique(n_z);

% n_z=L_el:L_el:H-L_el;
griglia_pareti=[];
for ij=1:numel(n_z);
   griglia_pareti_ij=griglia_bordo;
   griglia_pareti_ij(:,3)=n_z(ij)*ones(numel(griglia_bordo(:,1)),1);
   griglia_pareti=[griglia_pareti; griglia_pareti_ij];             
end
ndlist9=[[ndlist8(end,1)+1:ndlist8(end,1)+numel(griglia_pareti(:,1))]' griglia_pareti];


regular_mesh=[ndlist1; ndlist2];
circular_mesh=[ndlist2; ndlist3; ndlist4];
regular_mesh_top=[ndlist5; ndlist6];
circular_mesh_top=[ndlist6; ndlist7; ndlist8];
wall_mesh=[ndlist4; ndlist8; ndlist9];

for ij=1:numel(wall_mesh(:,1))
    alp=atan2(wall_mesh(ij,3),wall_mesh(ij,2));           
    H_ij=wall_mesh(ij,4);   
    mapped_wm(ij,:)=[wall_mesh(ij,1) alp H_ij];
end

mapped_wm_plusz=mapped_wm(mapped_wm(:,2)>=0,:);
mapped_wm_minz=mapped_wm(mapped_wm(:,2)<=0,:);

alpmax=max(mapped_wm_plusz(:,2));
alpmin=min(mapped_wm_minz(:,2));

mapped_wm_sel1=mapped_wm_plusz(mapped_wm_plusz(:,2)==alpmax,:);
mapped_wm_sel2=mapped_wm_minz(mapped_wm_minz(:,2)==alpmin,:);
mapped_wm_sel=[mapped_wm_sel1;mapped_wm_sel2];

% % selezione nodi per rinforzi in tubolare verticali
% for ij=1:numel(alpR)
%     alp=alpR(ij);
%     vert_sel=(abs(mapped_wm(:,2)/alp)<1.01).*(abs(mapped_wm(:,2)/alp)>0.99);
%     index_vert_sel=find(vert_sel==1);
%     
%     allbeam2{ij}=wall_mesh(index_vert_sel,:);
% end


% Fondo serbatoio - proprietà materiale 1
quad=quadronode(regular_mesh);
quad=[(1:numel(quad(:,1)))' ones(numel(quad(:,1)),1) ones(numel(quad(:,1)),1) quad];
polarquad=polarnode(circular_mesh);
polarquad=[(quad(end,1)+1:quad(end,1)+numel(polarquad(:,1)))'...
    2*ones(numel(polarquad(:,1)),1) ones(numel(polarquad(:,1)),1) polarquad];
quad(:,8)=zeros(numel(quad(:,1),1));
polarquad(:,8)=zeros(numel(polarquad(:,1),1));

% Top serbatoio - proprietà materiale 2 
quadtop=quadronode(regular_mesh_top);
quadtop=[(polarquad(end,1)+1:polarquad(end,1)+numel(quadtop(:,1)))' ...
    3*ones(numel(quadtop(:,1)),1) 2*ones(numel(quadtop(:,1)),1) quadtop];
polarquadtop=polarnode(circular_mesh_top);
polarquadtop=[(quadtop(end,1)+1:quadtop(end,1)+numel(polarquadtop(:,1)))' ...
    4*ones(numel(polarquadtop(:,1)),1) 2*ones(numel(polarquadtop(:,1)),1) polarquadtop];

% pareti serbatoio - proprietà materiale 3
wallquad=quadronode(mapped_wm_plusz);
wallquad=[(polarquadtop(end,1)+1:polarquadtop(end,1)+numel(wallquad(:,1)))' ...
    5*ones(numel(wallquad(:,1)),1) 3*ones(numel(wallquad(:,1)),1) wallquad];

wallquad_f=quadronode(mapped_wm_minz);
wallquad_f=wallquad_f(:,[4 3 2 1 5]);
wallquad_f=[(wallquad(end,1)+1:wallquad(end,1)+numel(wallquad_f(:,1)))' ...
    6*ones(numel(wallquad_f(:,1)),1) 4*ones(numel(wallquad_f(:,1)),1) wallquad_f];

wallquad_f2=quadronode(mapped_wm_sel);
wallquad_f2=wallquad_f2(:,[4 3 2 1 5]);
wallquad_f2=[(wallquad_f(end,1)+1:wallquad_f(end,1)+numel(wallquad_f2(:,1)))' ...
    7*ones(numel(wallquad_f2(:,1)),1) 4*ones(numel(wallquad_f2(:,1)),1) wallquad_f2];

ndlist=[ndlist1;ndlist2;ndlist3;ndlist4;ndlist5;ndlist6;ndlist7;ndlist8;ndlist9];
quad=[quad ;polarquad; quadtop ; polarquadtop; wallquad; wallquad_f; wallquad_f2];   

%% BEAM
% inspessimento serbatoio
iz=1;
for ij=1:numel(y_c);
    node4beam=wall_mesh(wall_mesh(:,4)==y_c(ij),:);
        
        for ik=1:numel(node4beam(:,1))   
            alp=atan2(node4beam(ik,3),node4beam(ik,2));
            if alp<0
               alp=alp+2*pi; 
            end
            B_list(ik,:)=[node4beam(ik,1) alp];
        end
    [unused1,isort]=sort(B_list(:,2));
    B_list=B_list(isort,:);
    clear beam_nodes
    for ik=1:numel(B_list(:,1))-1
        beam_nodes(ik,:)=[iz 1 1 B_list(ik,1) B_list(ik+1,1)];
        iz=iz+1;
    end        
        beam_nodes(ik+1,:)=[iz 1 1 B_list(end,1) B_list(1,1)];
        iz=iz+1;
        allbeam{ij}=beam_nodes;  
end

% rinforzo con tubolare
finalvect=allbeam{end};
iz=finalvect(end,1)+1;
ij_start=numel(allbeam)+1;
for ij=1:numel(y_cr);
    node4beam=wall_mesh(wall_mesh(:,4)==y_cr(ij),:);
        
        for ik=1:numel(node4beam(:,1))   
            alp=atan2(node4beam(ik,3),node4beam(ik,2));
            if alp<0
               alp=alp+2*pi; 
            end
            B_list(ik,:)=[node4beam(ik,1) alp];
        end
    [unused2,isort]=sort(B_list(:,2));
    B_list=B_list(isort,:);
    
     clear beam_nodes
    for ik=1:numel(B_list(:,1))-1
        %     id gruppo proprietà nodo 1 nodo 2   
        beam_nodes(ik,:)=[iz 1 2 B_list(ik,1) B_list(ik+1,1)];
        iz=iz+1;
    end        
        beam_nodes(ik+1,:)=[iz 1 2 B_list(end,1) B_list(1,1)];
        iz=iz+1;
        allbeam{ij_start}=beam_nodes;  
        ij_start=ij_start+1;
end

% selezione nodi per rinforzi in tubolare verticali
finalvect=allbeam{end};
iz=finalvect(end,1)+1;
ij_start=numel(allbeam)+1;
top_R=max(y_cr);


for ij=1:numel(alpR)
    
    alp=alpR(ij);
    if alp==0
    vert_sel=((mapped_wm(:,2))==0).*(mapped_wm(:,3)<=top_R);   
    else
    vert_sel=((mapped_wm(:,2)/alp)<1.01).*((mapped_wm(:,2)/alp)>0.99).*(mapped_wm(:,3)<=top_R);
    end
    
    index_vert_sel=find(vert_sel==1);
    keepindexv{ij}=index_vert_sel;
 
    node4beam=wall_mesh(index_vert_sel,:);
    B_list=node4beam;
%             for ik=1:numel(node4beam(:,1))   
%             alp=atan2(node4beam(ik,3),node4beam(ik,2));
%             if alp<0
%                alp=alp+2*pi; 
%             end
%             B_list(ik,:)=[node4beam(ik,1) alp];
%             end
            
                [unused2,isort]=sort(B_list(:,2));
                B_list=B_list(isort,:);
     clear beam_nodes
    if isempty(B_list)~=1
    for ik=1:numel(B_list(:,1))-1
        %     id gruppo proprietà nodo 1 nodo 2   
        beam_nodes(ik,:)=[iz 1 3 B_list(ik,1) B_list(ik+1,1)];
        iz=iz+1;
    end        
%         beam_nodes(ik+1,:)=[iz 1 3 B_list(end,1) B_list(1,1)];
%         iz=iz+1;
        allbeam{ij_start}=beam_nodes;  
        orient_vect(ij_start,:)=[ij_start alp];
        ij_start=ij_start+1;
    end
    
end


%% MODEL INFORMATION

    
    pn=fopen(nome_modello,'w');
    
% Scrittura dell prime informazioni generali sul file di testo, usate da
% Straus per realizzare poi il modello
fprintf(pn,'%s %s\n','FileFormat','Straus7.2.3.3');
fprintf(pn,'%s %s\n','ModelName','"Modello_serbatoio_cilindrico"');
fprintf(pn,'%s %s\n','Title','""');
fprintf(pn,'%s %s\n','Project','""');
fprintf(pn,'%s %s\n','Author','""');
fprintf(pn,'%s %s\n','Reference','""');

    
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

% Load Case 1
fprintf(pn,'%s %s %s %s\n','LoadCase','1','1','"Gravity"');
fprintf(pn,'%s %s %s\n', 'Gravity','3','-9.81000000000000E+0');
fprintf(pn,'%s %s\n','LCInclude','3');

% Load Case 2
fprintf(pn,'%s %s %s %s\n','LoadCase','2','0','"Serbatoio_pressione"');
fprintf(pn,'%s %s\n','LCInclude','3');

% Load Case 3
fprintf(pn,'%s %s %s %s\n','LoadCase','3','0','"Serbatoio_depressione"');
fprintf(pn,'%s %s\n','LCInclude','3');

% Load Case 4
fprintf(pn,'%s %s %s %s\n','LoadCase','4','0','"Vento"');
fprintf(pn,'%s %s\n','LCInclude','3');

% Load Case 5
fprintf(pn,'%s %s %s %s\n','LoadCase','5','0','"Uomo"');
fprintf(pn,'%s %s\n','LCInclude','3');


%% COORDINATE SYSTEM DEFINITIONS

% Definizione del sistema di coordinate principale
fprintf(pn,'%s %s %s %s\n','CoordSys','1','"Cartesiano"','GlobalXYZ');
fprintf(pn,'%s %s %s %s\n','CoordSys','2','"Polare"','PolarXYZ 0 0 0');



% % Stampa della matrice dei nodi su file txt
for i=1:numel(ndlist(:,1));
% Nell'ordine va posto l'id del nodo e le coordinate dello stesso in X,Y,Z
fprintf(pn,'%s %s %s %s %s\n','Node',num2str(ndlist(i,1)),...
     num2str(ndlist(i,2)),num2str(ndlist(i,3)),num2str(ndlist(i,4)));
end

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

for ij=1:numel(allbeam)
    if orient_vect(ij,1)~=0
        beam=allbeam{ij};
        angle=90+orient_vect(ij,2)/pi*180;
       
    for ik=1:numel(beam(:,1))        
        fprintf(pn,'%s %s %s\n','BmAngle',num2str(beam(ik,1)),num2str(angle));
        jk=jk+1;
    end
    end
end


for ij=1:numel(allbeam)
    beam=allbeam{ij};
        
        if beam(1,3)==1
        offset1=0;
        offset2=(-D/2-s_p/2);   
        else
        offset1=0;
        offset2=(-Dr/2-s_p/2);    
        end
    
    for ik=1:numel(beam(:,1))        
        fprintf(pn,'%s %s %s %s %s\n','BmOffset',num2str(beam(ik,1)),num2str(offset1),'  ',num2str(offset2));
        jk=jk+1;
    end
end

for i=1:numel(quad(:,1));
    % Nell'ordine va posto l'id del nodo, il gruppo a cui appartiene  
    % l'elemento, la proprietà e i 4 nodi che lo compongono
    if quad(i,2)==1 || quad(i,2)==4 || quad(i,2)==6 
    fprintf(pn,'%s %s %s %s %s %s %s %s\n','Quad4',...
            num2str(quad(i,1)),num2str(quad(i,2)),num2str(quad(i,3)),...
            num2str(quad(i,7)),num2str(quad(i,6)),num2str(quad(i,5)),num2str(quad(i,4)));    
    else
    fprintf(pn,'%s %s %s %s %s %s %s %s\n','Quad4',...
            num2str(quad(i,1)),num2str(quad(i,2)),num2str(quad(i,3)),...
            num2str(quad(i,4)),num2str(quad(i,5)),num2str(quad(i,6)),num2str(quad(i,7)));
    end
end

for i=1:numel(quad(:,1));
%     % Nell'ordine va posto l'id dell'elemento e il valore di pressione
%     % corrispondente
    if y_h2o-quad(i,8)<0 || quad(i,3)==2
    p=-pr_max;    
    else
    p=-9.806*density*(y_h2o-quad(i,8))*10^-6-pr_max;
    end
        fprintf(pn,'%s %s %s %s %s\n','Load Case 2','PlPressure',...
                '2',num2str(quad(i,1)),...
                num2str(p));
    
end

for i=1:numel(quad(:,1));
%     % Nell'ordine va posto l'id dell'elemento e il valore di pressione
%     % corrispondente
    
    p=-pr_min;    
    
        fprintf(pn,'%s %s %s %s %s\n','Load Case 3','PlPressure',...
                '3',num2str(quad(i,1)),...
                num2str(p));
    
end

for i=1:numel(quad(:,1));
%     % Nell'ordine va posto l'id dell'elemento e il valore di pressione
%     % corrispondente
    if quad(i,3)==3
        p=pwind_front;    
    
        fprintf(pn,'%s %s %s %s %s %s %s\n','Load Case 4','PlGlobalLoad',...
                '4',num2str(quad(i,1)),'0',...
                num2str(p),'0');        
    end
    if quad(i,3)==4
        p=pwind_back;    
    
        fprintf(pn,'%s %s %s %s %s %s %s\n','Load Case 4','PlGlobalLoad',...
                '4',num2str(quad(i,1)),'0',...
                num2str(p),'0');        
    end    
end





%% BEAM PROPERTIES

% Beam Property 1: inspessimento polipropilene
    A=B*D;
    I11=1/12*(B*D^3);
    I22=1/12*(B^3*D);
    J=2*1/3*B*D;


    fprintf(pn,'%s %s\n','BeamProp','1','"Spessore_PP"');
    fprintf(pn,'%s %s\n','MaterialName','"PP"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','UsePoisson','TRUE');
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','Expansion',num2str(10));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));
    fprintf(pn,'%s %s\n','Area',num2str(A));    
    fprintf(pn,'%s %s\n','MomentI11',num2str(I11));
    fprintf(pn,'%s %s\n','MomentI22',num2str(I22));
    fprintf(pn,'%s %s\n','MomentJ',num2str(J));
    fprintf(pn,'%s %s\n','SectionType','SolidRect');
    fprintf(pn,'%s %s\n','B',num2str(B));
    fprintf(pn,'%s %s\n','D',num2str(D));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');

    % Beam Property 2: Rinforzi di tubolare
    A=Br*Dr-(Br-2*T1r)*(Dr-2*T1r);
    I11=1/12*(Br*Dr^3)-1/12*((Br-2*T1r)*(Dr-2*T2r)^3);
    I22=1/12*(Br^3*Dr)-1/12*((Br-2*T1r)^3*(Dr-2*T2r));
    J=2*1/3*Br*T1r^3+2*1/3*Dr*T2r^3;


    fprintf(pn,'%s %s\n','BeamProp','2','"Cerchiature"');
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
    fprintf(pn,'%s %s\n','MomentJ',num2str(J));
    fprintf(pn,'%s %s\n','SectionType','HollowRect');
    fprintf(pn,'%s %s\n','B',num2str(Br));
    fprintf(pn,'%s %s\n','D',num2str(Dr));
    fprintf(pn,'%s %s\n','T1',num2str(T1r));
    fprintf(pn,'%s %s\n','T2',num2str(T2r));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    
        % Beam Property 3: Rinforzi di tubolare
    A=Br*Dr-(Br-2*T1r)*(Dr-2*T1r);
    I11=1/12*(Br*Dr^3)-1/12*((Br-2*T1r)*(Dr-2*T2r)^3);
    I22=1/12*(Br^3*Dr)-1/12*((Br-2*T1r)^3*(Dr-2*T2r));
    J=2*1/3*Br*T1r^3+2*1/3*Dr*T2r^3;


    fprintf(pn,'%s %s\n','BeamProp','3','"Rinforzi_vert"');
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
    fprintf(pn,'%s %s\n','MomentJ',num2str(J));
    fprintf(pn,'%s %s\n','SectionType','HollowRect');
    fprintf(pn,'%s %s\n','B',num2str(Br));
    fprintf(pn,'%s %s\n','D',num2str(Dr));
    fprintf(pn,'%s %s\n','T1',num2str(T1r));
    fprintf(pn,'%s %s\n','T2',num2str(T2r));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');

%% PLATE PROPERTIES

% Proprietà per gli elementi del fondo
    fprintf(pn,'%s %s %s\n','PlateShellProp','1','"Fondo"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_f));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_f));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');
    
    % Proprietà per gli elementi del top    
    fprintf(pn,'%s %s %s\n','PlateShellProp','2','"Top"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));    
    fprintf(pn,'%s %s\n','MemThick',num2str(s_t));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_t));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');   
    
% Proprietà per gli elementi delle pareti
    fprintf(pn,'%s %s %s\n','PlateShellProp','3','"Pareti1"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_p));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_p));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');    
 
% Proprietà per gli elementi delle pareti
    fprintf(pn,'%s %s %s\n','PlateShellProp','5','"Pareti2"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_p));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_p));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');    
    
    % Proprietà per gli elementi delle pareti
    fprintf(pn,'%s %s %s\n','PlateShellProp','3','"\\Pareti"');
    fprintf(pn,'%s %s\n','MaterialName','"mat"');
    fprintf(pn,'%s %s\n','Modulus',num2str(E_pp));
    fprintf(pn,'%s %s\n','Poisson',num2str(ni_pp));
    fprintf(pn,'%s %s\n','Density',num2str(rho_pp));
    fprintf(pn,'%s %s\n','ThermalCond',num2str(cond_pp));
    fprintf(pn,'%s %s\n','SpecificHeat',num2str(spec_heat_pp));
    fprintf(pn,'%s %s\n','MemThick',num2str(s_p));
    fprintf(pn,'%s %s\n','BendThick',num2str(s_p));
    fprintf(pn,'%s %s\n','NonLinType','Elasticplastic');
    fprintf(pn,'%s %s\n','YieldCriterion','VonMises');
    fprintf(pn,'%s %s\n','NumLayers','10');    

fclose all;


