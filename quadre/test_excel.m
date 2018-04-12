% lettura parametri vasche da file Excel

% per octave:
% pkg install C:/Users/sae/Downloads/io-2.4.10.tar.gz
% pkg load io

clear all
close all
pkg load io

[sheet_numbers,sheet_text,sheet_raw]=xlsread('Input_vasca.xlsx');

headers=sheet_text(1,:);

righe_lette=sheet_raw(:,1);

n_vasche=numel(sheet_raw(:,1));

count_nan=0;

for ik=1:n_vasche
    if isnan(righe_lette{ik})
        count_nan=count_nan+1;
    end
end

n_vasche=n_vasche-count_nan;

trova_indice(headers, 'Nome')
trova_indice(headers, 'path')