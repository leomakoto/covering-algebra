%Recebe: Os S premontados, [p;q], o grau dos S (d), o grau dos q (dq), o
%vetor de variaveis var e o vetor de centros simbolicos c.

%Retorna: O jacobiano de coeffs(Psatz), na seguinte ordem de derivadas: 
%1) coordenada dos centros (c11,c21,c31,...)
%2) indice dos centros (c11,c21,c31,c12,c22,c32,...) para n=3.
%3) raio
%4) triangular inferior das S com diagonal, de cima pra baixo, esquerda-direita.
%5) indices das S.

function J=jacobyalmip(S,pq,d,dq,var,c)
    m=monvecyalmip(var,d); M=m*m.';
    n=length(var);
    nc=length(c(1,:));
    ns=length(S(:,1));
    s=length(S(1,:))/ns;
    nq=length(dq);
    tr=sdpvar(s,1);
    for i=1:s
        tr(i)=trace(M*S(:,((i-1)*ns+1):(i*ns)));
    end
    tam=nc*n+1;
    J=sdpvar(1,tam);
    
    %Derivadas nos centros.
    for i=1:nc
        for j=1:n
            J(n*(i-1)+j)=tr(i+1)*2*(c(j,i)-var(j));
        end
    end
    
    %Derivada no raio.
    J(tam)=-[0,ones(1,s-nq-1),zeros(1,nq)]*tr;
    Js=sdpvar(ns*(ns+1)*s/2,1).';
    pq=[1;pq];
    
    %Derivadas nos S.
    for i=1:s
        for j=1:ns
            for k=1:j
                %Triangular inferior + diagonal, e->d, c->b.
                aux=pq(i)*M(j,k);
                if j~=k
                    aux=aux*2;
                end
                Js((i-1)*ns*(ns+1)/2+j*(j-1)/2+k)=aux;
            end
        end
    end
    J=[J,Js];

    %Até aqui, temos o gradiente de Psatz, mas não de cada um dos coeffs.
    %Pegando os coeffs do gradiente, temos os gradientes dos coeffs.
    
    dp=2*d+max([2;dq]);
    Jtotal=sdpvar(nchoosek(n+dp,dp),length(J));

    m=monvecyalmip(var,dp);
    nm=length(m);
    [~,m]=coefficients(ones(1,nm)*m, var);

    %Melhorar para aproveitar esparsidade, talvez fazer uma função que,
    %dado os graus, retorna a posição em m, igual na deglex().
    for i=1:length(J)
        [auxc,auxv]=coefficients(J(i),var);
        for j=1:nm
            %Para depois q acabar o auxv.
            if length(auxv)<j
                auxv=vertcat(auxv,m(j));
                auxc=vertcat(auxc,0);
            %Para encaixar o zero onde nao tem, pois coefficients() retorna uma coisa sem zeros.
            elseif ~isequal(auxv(j),m(j))
                auxv=vertcat(auxv(1:(j-1)),m(j),auxv(j:end));
                auxc=vertcat(auxc(1:(j-1)),0,auxc(j:end));
            end
        end
        Jtotal(:,i)=auxc;
    end
    
    %Retorna o Jacobiano total dos coeffs.
    J=Jtotal;
end