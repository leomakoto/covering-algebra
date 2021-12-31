%Recebe: A matriz de restrições Cons, os coeffs do psatz S, o vetor de centros simbolicos c, o
%raio simbolico r, os centros e raio numericos numc, um raio inicial.

%Retorna: um raio inviável e um viável.

function [rf,ri]=factyalmip(Cons,S,c,r,numc,rini)
    [stat,~]=iscoveryalmip(Cons,S,c,r,numc,rini);
    if stat.problem~=1
        rf=rini;
        ri=0;
    else
        k=0;
        while stat.problem==1
            k=k+1;
            [stat,~]=iscoveryalmip(Cons,S,c,r,numc,rini+k*2^k);
        end        
        rf=rini+k*2^k;
        ri=rini+(k-1)*2^(k-1);
    end
end