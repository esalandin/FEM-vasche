function out=nodes(x,Dx)
% --------------------------------------------------
% ultima modifica: 17/10/2017
% --------------------------------------------------
    x=unique(x);
    out1=[];
for ij=1:numel(x)-1
    out1=[out1 x(ij):(x(ij+1)-x(ij))/round(((x(ij+1)-x(ij)))/Dx):x(ij+1)];    
end
    out=unique(out1);
end
