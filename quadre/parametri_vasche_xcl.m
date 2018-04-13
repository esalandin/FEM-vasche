% estraggo da una riga di valori i parametri di una vasca

function par_vasca= parametri_vasche_xcl(riga_nomi,riga_valori)

% l'operatore {} serve per estrarre il contenuto della cella (come {1,1} per una cella 1x1)
par_vasca.nome= riga_valori(trova_indice(riga_nomi,'Nome')) {};
par_vasca.path= riga_valori(trova_indice(riga_nomi,'path')) {};
par_vasca.descr= riga_valori(trova_indice(riga_nomi,'descrizione')) {};
par_vasca.L_x= riga_valori(trova_indice(riga_nomi,'X_vasca [m]')) {};
par_vasca.L_y= riga_valori(trova_indice(riga_nomi,'Y_vasca [m]')) {};
par_vasca.L_z= riga_valori(trova_indice(riga_nomi,'Z_vasca [m]')) {};
par_vasca.s_f= riga_valori(trova_indice(riga_nomi,'Spessore fondo [m]')) {};
par_vasca.s_p= riga_valori(trova_indice(riga_nomi,'Spessore Pareti [m]')) {};
% per sb viene utilizzato lo stesso valore di s_p
par_vasca.sb= par_vasca.s_p;
par_vasca.B= riga_valori(trova_indice(riga_nomi,'dimensione maggiore tubolare di rinforzo [m]')) {};
par_vasca.D= riga_valori(trova_indice(riga_nomi,'dimensione minore tubolare di rinforzo [m]')) {};
par_vasca.tt= riga_valori(trova_indice(riga_nomi,'Spessore tubolare di rinforzo [m]')) {};
par_vasca.rho_l= riga_valori(trova_indice(riga_nomi,'densita liquido [kg/m3]')) {};

switch char(riga_valori(trova_indice(riga_nomi,'Materiale')))
    case 'Acciaio'
        out{12}=2;    
        par_vasca.mat=2;
    case 'PP'
        out{12}=1;   
       par_vasca.mat=1; 
    case 'Titanio'
        out{12}=3;
        par_vasca.mat=3;
    case 'PE'
        out{12}=4;
        par_vasca.mat=4;
end

par_vasca.c_x= abate(riga_valori(trova_indice(riga_nomi,'costole in x [m]')));
par_vasca.c_y= abate(riga_valori(trova_indice(riga_nomi,'cerchiatura in y [m]')));
par_vasca.c_z= abate(riga_valori(trova_indice(riga_nomi,'costole in z [m]')));
par_vasca.ad_nd_x= abate(riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in x [m]')));
par_vasca.ad_nd_z= abate(riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in z [m]')));

par_vasca.dx= riga_valori(trova_indice(riga_nomi,'dimensione target mesh x [m]')) {};
par_vasca.dy= riga_valori(trova_indice(riga_nomi,'dimensione target mesh y [m]')) {};
par_vasca.dz= riga_valori(trova_indice(riga_nomi,'dimensione target mesh z [m]')) {};

end

% trasforma una stringa con numeri separati da trattini in un'array
% esempio '1-2.3 3,4' viene convertito in [1.0000   2.300   3.4000]
function v=abate(c)
    if iscell(c)
        content= c{1};
    else
        content= c;
    endif
    if ischar(content)
      str= strrep(content, ',', '.'); % sostituiamo le virgole con punto per non avere problemi di locale
      str= strrep(str, '-', ' '); % sostituiamo - con spazio come separatore
      [v, ok_conversione]= str2num(str);
      if ~ok_conversione
        disp(['conversione "' c '" fallita']);
      endif
    else
      v=content;
    endif
endfunction