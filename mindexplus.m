%Recebe: o numero de variaveis n e o grau total d.
%Retorna: o vetor de indices de potencias em deglex, cuja coluna eh o indice da variavel e a linha eh o monomio.

function ind=mindexplus(n,d)
    aux=unique(nchoosek(repmat(0:d,1,n),n),'rows');
    ind=deglex(aux,d);
end