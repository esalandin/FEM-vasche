function indice= trova_indice(celle, nome_campo)
  find_ret= find(strcmp(celle, nome_campo));
  if isempty(find_ret)
    indice=0;
   else
    indice= find_ret;    
  endif
endfunction
