%Recebe:
% linCons = psatz linearizado simbolico
% d = direcao simbolica
% c = centros simbolicos
% r = raio simbolico
% S = psatz-vars simbolicas
% numc = centros numericos viaveis
% numr = raio numerico viavel
% numS = psatz-vars numericas viaveis (da fase de restauracao)
% locr = posicao do raio em d

%Retorna:
% stat = resultado da otimizacao.
% dir = direcao de descida no espaco tangente a Cons saindo de (numc,r,numS).
% Se dir=0, o ponto atual é otimo (será?).

function [dir, stat]=tanDirYalmip(linCons, d, c, r, S, D, dcmat, numc, numr, numS, locr, sigma)
    linCons = replace(linCons, c, numc);
    linCons = replace(linCons, r, numr);
    linCons = replace(linCons, S, numS);
    %Duvidoso, conserte.
    
    %norm(D,inf)<=sigma*norm(numS,inf); norm(dcmat,inf)<=sigma*norm(numc,inf); 
    stat=optimize(linCons+[d(locr)+numr>=0;norm(D,inf)<=sigma*norm(numS,inf)],numr+d(locr),sdpsettings('verbose',0,'solver','sedumi','sedumi.eps',1e-12));
    dir = value(d);
end