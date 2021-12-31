%Obs: complemento ao mindexplus.
%Recebe: uma matriz de indices que em cada linha tem um vetor de potencias e o grau max.
%Retorna: o vetor ordenado em deglex e remove os de grau improprio > d.

function u=deglex(v,d)
    %n=length(v(:,1));
    m=length(v(1,:));
    on=ones(1,m).*(100^(m+1));
    cr=100.^(1:m);
    %deg de cada linha.
    deg=(on*v.').';
    %peso de cada potencia.
    pes=(cr*v.').';
    [~,I]=sort(deg+pes);
    v=v(I,:);
    fnd=find(sort(deg)==d*100^(m+1));
    u=v(1:fnd(end),:);
end