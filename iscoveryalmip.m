%Recebe:
% Cons = matriz de restricoes do psatz
% S = coeffs do psatz
% c = vet de centros simbolicos
% r = raio simbolico
% numc = centro numerico
% numr = raio numerico

%Retorna:
% stat = se c=numc, r=numr forma um SDP viável ou não
% vS = o valor de S

function [stat,vS]=iscoveryalmip(Cons,S,c,r,numc,numr)
    Cons=replace(Cons,c,numc);
    Cons=replace(Cons,r,numr);
    stat=optimize(Cons,norm(S,inf),sdpsettings('verbose',0,'solver','sedumi'));
    vS=value(S);
end