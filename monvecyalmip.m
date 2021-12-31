%Recebe: vetor de variaveis ordenado x1,x2,...,xn (var), grau maximo dos
%monomios d.

%Retorna: o vetor de monomios de grau ate d ordenado em deglex, x1<...<xn.

%Observação: A funcao coefficients() do Yalmip usa lex com xn<...<x1. Para
%colocar m na mesma ordem, faca m=coefficients(ones.'*m,var);
function m=monvecyalmip(var,d)
    n=length(var);
    I=mindexplus(n,d);
    k=nchoosek(n+d,d);
    m=sdpvar(k,1);
    for i=1:k
        m(i)=prodyalmip(var.^(I(i,:).'));
    end
end