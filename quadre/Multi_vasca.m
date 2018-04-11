
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------

clear all
close all
clc

% global keepvalue keepvalue2

[n,t,r]=xlsread('Input_vasca_CR.xlsx');
headers=r(1,:);
righe_lette=r(:,1);

n_vasche=numel(r(:,1));

count_nan=0;

for ik=1:n_vasche
    if isnan(righe_lette{ik})
        count_nan=count_nan+1;
    end
end

n_vasche=n_vasche-count_nan;

for ij=1:n_vasche-1
    
%     nome          par{1}      identificativo vasca 
%     L_x           par{4}      lunghezza vasca
%     L_y           par{5}      profondità vasca
%     L_z           par{6}
%     s_f           par{7}
%     s_p           par{8}
%     sb            default=spessore pareti
%     B             par{9}
%     D             par{10}        
%     t             par{11}
%     c_x           par{13}
%     c_y           par{14}    
%     c_z           par{15}
%     mat           par{12}
%     descrizione   par{2}
%     path          par{3}
%     rho_l         par{16}  
%     apg           par{17}
%     pos_apg       par{18}   
%     ad_nd_x       par{19}             nodi di interesse in aggiunta su asse x    
%     ad_nd_z       par{20}             nodi di interesse in aggiunta su asse z    

    

    par=parametri_vasche_xcl(headers(1,:),r(ij+1,:));
%     par{17}='HEA100';
%     par{18}=[0.06 1.16];
%     par{19}=
%     par{20}=
    
    VG(par{1},par{4},par{5},par{6},par{7},par{8},par{8},par{9},par{10},par{11},par{13},par{14},par{15},par{12},par{2},par{3},par{16},par{17},par{18},par{19},par{20})
end