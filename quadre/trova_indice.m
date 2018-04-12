function indice= trova_indice(celle, nome_campo)
  nome_campo
  find_ret= find(strcmp(celle, nome_campo));
  if isempty(find_ret)
    disp(['trova_indice ' '"' nome_campo '"' ' non trovato']);
   else
    disp(['trova_indice ' '"' nome_campo '"' ' indice= ' num2str(find_ret)]);
  endif    
  if isempty(find_ret)
    indice=0
   else
    indice= find_ret    
  endif
endfunction
