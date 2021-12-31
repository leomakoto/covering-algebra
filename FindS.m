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

function [stat,vSr]=FindS(PosdefCons,S,H,Hp,r,Rp,par)
    Hpp=replace(H,r,Rp);
    %norm(S,inf)
    stat=optimize(PosdefCons+[norm(Hpp,inf)<=par*Hp],0,sdpsettings('verbose',0,'solver','sedumi'));
    vSr=value(S);
end