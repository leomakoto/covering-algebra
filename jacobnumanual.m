function numJ=jacobnumanual(H,S,c,r,nums,numc,numr,locr,delta)
    %Colunas e linhas de S.
    cs=length(S(1,:)); ls=length(S(:,1));
    %Numero de coeffs de S.
    ns=cs/ls;
    %Colunas de d para cada s.
    %cd=ls*(ls+1)/2;
    %Linhas e colunas de c.
    n=length(c(:,1)); nc=length(c(1,:));
    %Calcular a direção em c, r e S.
    tam=(nc*n)+length(r)+ls*(ls+1)*ns/2;
    d=sdpvar(1,tam);
    numJ=zeros(length(H),tam);
    for i=1:nc
        for j=1:n
            d(n*(i-1)+j)=c(j,i);
        end
    end
    %d((locr+1):(locr+nc))=r;
    d(locr)=r;
    for i=1:ns
        for j=1:ls
            for k=1:j
                %Triangular inferior + diagonal, e->d, c->b.
                d(locr+(i-1)*ls*(ls+1)/2+j*(j-1)/2+k)=S(j,ls*(i-1)+k);
            end
        end
    end
    numd=avalcrs(d,c,r,S,numc,numr,nums);
    numH=replace(H,d,numd);
    for i=1:length(d)
        numd(i)=numd(i)+delta;
        numJ(:,i)=(replace(H,d,numd)-numH)./delta;
        numd(i)=numd(i)-delta;
    end
end