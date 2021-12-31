%Recebe: A matriz de restrições Cons, os coeffs do psatz S, o vetor de centros simbolicos c, o
%raio simbolico r, os centros e raio numericos numc, rf e ri (viavel e
%inviavel) e precisão eps.

%Retorna: o raio ótimo para c=numc.

function [R,vS]=binaryalmip(Cons,S,c,r,numc,eps,rf,ri)
    dif=abs(rf-ri);
    R=rf;
    if dif>eps
        med=(rf+ri)/2;
        [stat,~]=iscoveryalmip(Cons,S,c,r,numc,med);
        numstat=stat.problem;
        if numstat==1
            [Raux,~]=binaryalmip(Cons,S,c,r,numc,eps,rf,med);
        else
            [Raux,~]=binaryalmip(Cons,S,c,r,numc,eps,med,ri);
        end
        R=Raux;
    end
    vS=value(S);
end