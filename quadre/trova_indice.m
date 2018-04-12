function indice= trova_indice(celle, nome_campo)
  disp('trova_indice');
  nome_campo
  indice= find(strcmp(celle, nome_campo));
  if isempty(indice)
    disp(['"' nome_campo '"' ' non trovato']);
   else
    disp(['"' nome_campo '"' ' indice= ' num2str(indice)]);
  endif    
endfunction
