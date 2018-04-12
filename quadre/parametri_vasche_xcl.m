function out=parametri_vasche_xcl(v1,v2)

out(1)=v2(trova_indice(v1,'Nome'));
out(2)=v2(trova_indice(v1,'path'));
out(3)=v2(trova_indice(v1,'descrizione'));
out(4)=v2(trova_indice(v1,'X_vasca [m]'));

out(5)=v2(trova_indice(v1,'Y_vasca [m]'));
out(6)=v2(trova_indice(v1,'Z_vasca [m]'));
out(7)=v2(trova_indice(v1,'Spessore fondo [m]'));
out(8)=v2(trova_indice(v1,'Spessore Pareti [m]'));
out(9)=v2(trova_indice(v1,'dimensione maggiore tubolare di rinforzo [m]'));
out(10)=v2(trova_indice(v1,'dimensione minore tubolare di rinforzo [m]'));
out(11)=v2(trova_indice(v1,'Spessore tubolare di rinforzo [m]'));
out(16)=v2(trova_indice(v1,'densita liquido [kg/m3]'));
temp_out(12)=v2(trova_indice(v1,'Materiale'));
    char12=char(temp_out{12});

switch char12
    case 'Acciaio'
        out{12}=2;    
    case 'PP'
        out{12}=1;    
    case 'Titanio'
        out{12}=3;
    case 'PE'
        out{12}=4;
end


temp_out(13)=v2(trova_indice(v1,'costole in x [m]'));
out{13}=abate(temp_out(13));     
    
temp_out(14)=v2(trova_indice(v1,'cerchiatura in y [m]'));
out{14}=abate(temp_out(14));  

%     char14=char(temp_out{14});
%     indexes=[0 regexp(char14,'-')];    
%     for ij=1:numel(indexes)-1;
%         v14(ij)=str2num(char14(indexes(ij)+1:indexes(ij+1)-1));
%     end
%     v14(ij+1)=str2num(char14(indexes(end)+1:end));
%     out{14}=sort(v14);        
    
temp_out(15)=v2(trova_indice(v1,'costole in z [m]'));
out{15}=abate(temp_out(15)); 

temp_out(19)=v2(trova_indice(v1,'posizioni aggiuntive in x [m]'));
out{19}=abate(temp_out(15)); 

temp_out(20)=v2(trova_indice(v1,'posizioni aggiuntive in z [m]'));
out{20}=abate(temp_out(20)); 
%     char15=char(temp_out{15});
%     indexes=[0 regexp(char15,'-')];
%     
%     for ij=1:numel(indexes)-1;
%         v15(ij)=str2num(char15(indexes(ij)+1:indexes(ij+1)-1));
%     end
%     v15(ij+1)=str2num(char15(indexes(end)+1:end));
%     out{15}=sort(v15);    
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
