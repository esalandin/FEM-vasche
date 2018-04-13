% lettura parametri vasche da file Excel

% per octave:
% pkg install C:/Users/sae/Downloads/io-2.4.10.tar.gz
% pkg load io

clear
close all
pkg load io

% global keepvalue keepvalue2

[n,t,r]=xlsread('Input_vasca_Z.xlsx');
headers=r(1,:);
righe_lette=r(:,1);

n_vasche=numel(r(:,1));

count_nan=0;

for ik=1:n_vasche
    if isnan(righe_lette{ik})
        count_nan=count_nan+1;
    endif
endfor

n_vasche=n_vasche-count_nan;

for ij=1:n_vasche-1
% par e' una struttura che contiene i parametri letti da una riga del file excel
    par=parametri_vasche_xcl(headers(1,:),r(ij+1,:));

% apg_type e apg_pos non vengono letti dal file excel
% valori di esempio:
% apg_type= 'HEA100';
% apg_pos= [0.06 1.16];
  par.apg_type= NaN;
  par.apg_pos= NaN;

  VG(par);
endfor
