function out= parametri_vasche_xcl(v1,v2)

% par{1},par{4},par{5},par{6},par{7},par{8},par{8},par{9},par{10},par{11},par{13},par{14},par{15},
% nome,  L_x,   L_y,   L_z,   s_f,   s_p,   sb,    B,     D,      tt,     c_x,    c_y,    c_z,    
% par{12},par{2},par{3},par{16},par{17} ,par{18},par{19},par{20}
% mat,    path,  descr, rho_l,  apg_type,apg_pos,ad_nd_x,ad_nd_z

% i parametri 17 e 18 (apg_type,apg_pos) non vengono letti dal fogliio excel.

out(1)=v2(trova_indice(v1,'Nome'));
par_struct.nome= v2(trova_indice(v1,'Nome'));

out(2)=v2(trova_indice(v1,'path'));
par_struct.path= v2(trova_indice(v1,'path'));

out(3)=v2(trova_indice(v1,'descrizione'));
par_struct.descr= v2(trova_indice(v1,'descrizione'));

out(4)=v2(trova_indice(v1,'X_vasca [m]'));
par_struct.L_x= v2(trova_indice(v1,'X_vasca [m]'));

out(5)=v2(trova_indice(v1,'Y_vasca [m]'));
par_struct.L_y= v2(trova_indice(v1,'Y_vasca [m]'));

out(6)=v2(trova_indice(v1,'Z_vasca [m]'));
par_struct.L_z= v2(trova_indice(v1,'Z_vasca [m]'));

out(7)=v2(trova_indice(v1,'Spessore fondo [m]'));
par_struct.s_f= v2(trova_indice(v1,'Spessore fondo [m]'));

out(8)=v2(trova_indice(v1,'Spessore Pareti [m]'));
par_struct.s_p= v2(trova_indice(v1,'Spessore Pareti [m]'));

% per sb viene utilizzato lo stesso valore di s_p
par_struct.sb= par_struct.s_p;

out(9)=v2(trova_indice(v1,'dimensione maggiore tubolare di rinforzo [m]'));
par_struct.B= v2(trova_indice(v1,'dimensione maggiore tubolare di rinforzo [m]'));

out(10)=v2(trova_indice(v1,'dimensione minore tubolare di rinforzo [m]'));
par_struct.D= v2(trova_indice(v1,'dimensione minore tubolare di rinforzo [m]'));

out(11)=v2(trova_indice(v1,'Spessore tubolare di rinforzo [m]'));
par_struct.tt= v2(trova_indice(v1,'Spessore tubolare di rinforzo [m]'));

out(16)=v2(trova_indice(v1,'densita liquido [kg/m3]'));
par_struct.rho_l= v2(trova_indice(v1,'densita liquido [kg/m3]'));

switch char(v2(trova_indice(v1,'Materiale')))
    case 'Acciaio'
        out{12}=2;    
        par_struct.mat=2;
    case 'PP'
        out{12}=1;   
       par_struct.mat=1; 
    case 'Titanio'
        out{12}=3;
        par_struct.mat=3;
    case 'PE'
        out{12}=4;
        par_struct.mat=4;
end


temp_out(13)=v2(trova_indice(v1,'costole in x [m]'));
out{13}=abate(temp_out(13));     
par_struct.c_x= abate(v2(trova_indice(v1,'costole in x [m]')));
    
temp_out(14)=v2(trova_indice(v1,'cerchiatura in y [m]'));
out{14}=abate(temp_out(14));  
par_struct.c_y= abate(v2(trova_indice(v1,'cerchiatura in y [m]')));

temp_out(15)=v2(trova_indice(v1,'costole in z [m]'));
out{15}=abate(temp_out(15)); 
par_struct.c_z= abate(v2(trova_indice(v1,'costole in z [m]')));

temp_out(19)=v2(trova_indice(v1,'posizioni aggiuntive in x [m]'));
out{19}=abate(temp_out(15)); 
par_struct.ad_nd_x= abate(v2(trova_indice(v1,'posizioni aggiuntive in x [m]')));

temp_out(20)=v2(trova_indice(v1,'posizioni aggiuntive in z [m]'));
out{20}=abate(temp_out(20)); 
par_struct.ad_nd_z= abate(v2(trova_indice(v1,'posizioni aggiuntive in z [m]')));

end

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
