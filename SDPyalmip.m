%Recebe:
% var = vetor de variaveis
% e = numero de esferas
% d = grau dos S
% q = vetor de polinomios do semialgebrico
% deg = grau do elemento diferente
% c = centros simbolicos
% r = raio simbolico

%Retorna:
% V = coeffs do psatz de Schmudgen
% b = (-1,0,0,...,0)
% S = psatz-vars simbolicos
% Cons = matriz de restricao montada
% circ = vetor de polinomios que def as esf

function [H,S,Cons,circ]=SDPyalmip(var,e,d,q,deg,c,r)
    n=length(var);
    m=length(q);
    p=sdpvar(e,1);
    for i=1:e
        p(i)=(var-c(:,i)).'*(var-c(:,i))-r;
    end
    ds=nchoosek(n+d,d);
    k=2^(e+m);
    Psatz=0;
    exp=dec2bin(0:(k-1));
    polvet=[p;q];
    v=monvecyalmip(var,d);
    S=sdpvar(ds,ds*k);
    s=sdpvar(k,1);
    for i=1:k
        S(:,((i-1)*ds+1):(i*ds))=sdpvar(ds,ds,'symmetric','real');
        s(i)=v.'*S(:,((i-1)*ds+1):(i*ds))*v;
        Psatz=Psatz+termopreorder(polvet,s(i),exp(i,:));
    end
    %DegPsatz=2*(e+d)+dq.'*ones(m,1);
    b=1;
    for i=1:e
       b=b*p(i)^deg;
    end
    Psatz=Psatz+b;
    H=coefficients(Psatz,var);
    Cons=[];
    for i=1:k
        Cons=Cons+[S(:,((i-1)*ds+1):(i*ds))>=zeros(ds,ds)];
    end
    circ=p;
end