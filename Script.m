%Numero de vars, esferas, nivel da hierarquia, grau dos qs.
n=2; e=1; d=2; dq=2;
symm=0;

%Centros e raios iniciais e precisão.
eps=10^(-5);
% numc=[1;1];
numc=rand(n,e);
%numc=[-0.5,0,0.5,-0.5,0,0.5;-0.5,-0.5,-0.5,0.5,0.5,0.5];
if symm==1
    numc(1,1)=0;
end

%Declaracao de vars.
var=sdpvar(n,1); sdpvar r; c=sdpvar(n,e,'full');

%Semialgebrico q>=0.
q=[-var(1)^2-var(2)^2+1.0];
%q=-var(1)^6-3*var(1)^4*var(2)^2+3*var(1)^4-3*var(1)^2*var(2)^4+var(1)^2*var(2)^3+6*var(1)^2*var(2)^2-3*var(1)^2-var(2)^6+3*var(2)^4-3*var(2)^2+1;
%const=sqrt(3/4);
%q=[var(2); const-var(2)-2*const*var(1); 2*const*var(1)-var(2)+const];

%Psatz de Putinar one-shot.
tic
%[H,S,PosdefCons,circ]=SDPyalmip(var,e,d,q,2,c,r);
[H,S,PosdefCons,circ]=SDPutyalmip(var,e,d,q,c,r);
pq=[circ;q];
Cons=PosdefCons+[H==zeros(length(H),1)];
toc

nr=[];
ns=[];
for i=1:10
[stat,vS]=iscoveryalmip(Cons,S,c,r,[0;0],1+1/i);
nr=[nr;1+1/(i*1000)];
ns=[ns;norm(vS)]; disp(i)
end

%Jacobiano one-shot simbolico.
% tic
% J=jacobyalmip(S,pq,d,dq,var,c);
% toc

%Solver.
tic
locr=e*n+1;
[rf,ri]=factyalmip(Cons,S,c,r,numc,0);
[Rp,vS]=binaryalmip(Cons,S,c,r,numc,10^(-2),rf,ri);
%[stat,vS]=iscoveryalmip(Cons,S,c,r,numc,Rp);
toc

s=0.001;
mdata=10;
dalpha=0;
sigma=1;
%tic; [Ropt,copt,I,M]=uniraioinexrest(PosdefCons, H, S, J, locr, c, r, numc, Rp, vS, eps, s, dalpha, symm, mdata, sigma); toc
tic; [Ropt,copt,I,M]=uniraioinexrestnum(PosdefCons, H, S, locr, c, r, numc, Rp, vS, eps, s, symm, mdata, sigma); toc
