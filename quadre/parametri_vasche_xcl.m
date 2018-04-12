function par_vasca= parametri_vasche_xcl(riga_nomi,riga_valori)

% par{1},par{4},par{5},par{6},par{7},par{8},par{8},par{9},par{10},par{11},par{13},par{14},par{15},
% nome,  L_x,   L_y,   L_z,   s_f,   s_p,   sb,    B,     D,      tt,     c_x,    c_y,    c_z,    
% par{12},par{2},par{3},par{16},par{17} ,par{18},par{19},par{20}
% mat,    path,  descr, rho_l,  apg_type,apg_pos,ad_nd_x,ad_nd_z

% i parametri 17 e 18 (apg_type,apg_pos) non vengono letti dal fogliio excel.
% l'operatore {} serve per estrarre il contenuto della cella (come {1,1} per una cella 1x1)

out(1)=riga_valori(trova_indice(riga_nomi,'Nome'));
par_vasca.nome= riga_valori(trova_indice(riga_nomi,'Nome')) {};

out(2)=riga_valori(trova_indice(riga_nomi,'path'));
par_vasca.path= riga_valori(trova_indice(riga_nomi,'path')) {};

out(3)=riga_valori(trova_indice(riga_nomi,'descrizione'));
par_vasca.descr= riga_valori(trova_indice(riga_nomi,'descrizione')) {};

out(4)=riga_valori(trova_indice(riga_nomi,'X_vasca [m]'));
par_vasca.L_x= riga_valori(trova_indice(riga_nomi,'X_vasca [m]')) {};

out(5)=riga_valori(trova_indice(riga_nomi,'Y_vasca [m]'));
par_vasca.L_y= riga_valori(trova_indice(riga_nomi,'Y_vasca [m]')) {};

out(6)=riga_valori(trova_indice(riga_nomi,'Z_vasca [m]'));
par_vasca.L_z= riga_valori(trova_indice(riga_nomi,'Z_vasca [m]')) {};

out(7)=riga_valori(trova_indice(riga_nomi,'Spessore fondo [m]'));
par_vasca.s_f= riga_valori(trova_indice(riga_nomi,'Spessore fondo [m]')) {};

out(8)=riga_valori(trova_indice(riga_nomi,'Spessore Pareti [m]'));
par_vasca.s_p= riga_valori(trova_indice(riga_nomi,'Spessore Pareti [m]')) {};

% per sb viene utilizzato lo stesso valore di s_p
par_vasca.sb= par_vasca.s_p;

out(9)=riga_valori(trova_indice(riga_nomi,'dimensione maggiore tubolare di rinforzo [m]'));
par_vasca.B= riga_valori(trova_indice(riga_nomi,'dimensione maggiore tubolare di rinforzo [m]')) {};

out(10)=riga_valori(trova_indice(riga_nomi,'dimensione minore tubolare di rinforzo [m]'));
par_vasca.D= riga_valori(trova_indice(riga_nomi,'dimensione minore tubolare di rinforzo [m]')) {};

out(11)=riga_valori(trova_indice(riga_nomi,'Spessore tubolare di rinforzo [m]'));
par_vasca.tt= riga_valori(trova_indice(riga_nomi,'Spessore tubolare di rinforzo [m]')) {};

out(16)=riga_valori(trova_indice(riga_nomi,'densita liquido [kg/m3]'));
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


temp_out(13)=riga_valori(trova_indice(riga_nomi,'costole in x [m]'));
out{13}=abate(temp_out(13));     
par_vasca.c_x= abate(riga_valori(trova_indice(riga_nomi,'costole in x [m]')));
    
temp_out(14)=riga_valori(trova_indice(riga_nomi,'cerchiatura in y [m]'));
out{14}=abate(temp_out(14));  
par_vasca.c_y= abate(riga_valori(trova_indice(riga_nomi,'cerchiatura in y [m]')));

temp_out(15)=riga_valori(trova_indice(riga_nomi,'costole in z [m]'));
out{15}=abate(temp_out(15)); 
par_vasca.c_z= abate(riga_valori(trova_indice(riga_nomi,'costole in z [m]')));

temp_out(19)=riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in x [m]'));
out{19}=abate(temp_out(15)); 
par_vasca.ad_nd_x= abate(riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in x [m]')));

temp_out(20)=riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in z [m]'));
out{20}=abate(temp_out(20)); 
par_vasca.ad_nd_z= abate(riga_valori(trova_indice(riga_nomi,'posizioni aggiuntive in z [m]')));

end

% trasforma una stringa con numeri separati da trattini in un'array
function v=abate(c)

    if iscell(c)
        c1=c{1};
    end
        if ischar(c1)
            ch=char(c1);
            indexes=[0 regexp(ch,'-')];          
            for ij=1:numel(indexes)-1;
                v(ij)=str2num(ch(indexes(ij)+1:indexes(ij+1)-1));
            end
            v(ij+1)=str2num(ch(indexes(end)+1:end));
            if indexes==0
                v_temp=char(c1);
                v=str2num(v_temp);
            end
        else
            v=c1;
        end
end
