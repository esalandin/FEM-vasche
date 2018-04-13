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
