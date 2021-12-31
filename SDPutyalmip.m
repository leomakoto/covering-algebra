%Recebe:
% var = vetor de variaveis
% e = numero de esferas
% d = grau dos S
% q = vetor de polinomios do semialgebrico
% dq = vetor de grau dos qs
% c = centros simbolicos
% r = raio simbolico

%Retorna:
% V = coeffs do psatz de Putinar
% b = (-1,0,0,...,0)
% S = psatz-vars simbolicos
% Cons = matriz de restricao montada
% circ = vetor de polinomios que def as esf

%funcoes boas: sdisplay, assign, value, replace.

function [H,S,Cons,circ]=SDPutyalmip(var,e,d,q,c,r)
    n=length(var);
    m=length(q);
    p=sdpvar(e,1);
    for i=1:e
        p(i)=(var-c(:,i)).'*(var-c(:,i))-r;
    end
    ds=nchoosek(n+d,d);
    k=1+e+m;
    Psatz=1;
%     Psatz=0;
    polvet=[1;p;q];
    v=monvecyalmip(var,d);
    S=sdpvar(ds,ds*k);
    s=sdpvar(k,1);
    for i=1:k
        S(:,((i-1)*ds+1):(i*ds))=sdpvar(ds,ds,'symmetric','real');
        s(i)=v.'*S(:,((i-1)*ds+1):(i*ds))*v;
        Psatz=Psatz+polvet(i)*s(i);
    end
%     DegPsatz=2*d+max([2;dq]);
%     V=coefficients(Psatz,var);
%     b=zeros(nchoosek(n+DegPsatz, DegPsatz),1);
%     b(1)=-1;
    H=coefficients(Psatz, var);
    Cons=[];
    for i=1:k
        Cons=Cons+[S(:,((i-1)*ds+1):(i*ds))>=zeros(ds,ds)];
    end
    circ=p;
end