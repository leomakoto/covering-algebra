function [Ropt,copt,I]=uniraioinexrest(PosdefCons, V, b, S, J, locr, c, r, numc, Rp, vS, eps, sigma, dalpha)           
    %Colunas de S.
    cs=length(S(1,:));
    %Linhas de S.
    ls=length(S(:,1));
    %Numero de coeffs de S.
    ns=cs/ls;
    %Linhas de c.
    lc=length(c(:,1));
    %Colunas de c.
    cc=length(c(1,:));
    %Matrizona D;
    D=sdpvar(ls,cs);
    
    %Direção nos S.
    for j=1:ns
        D(:,((j-1)*ls+1):(j*ls))=sdpvar(ls,ls,'symmetric','real');
    end
    d=sdpvar(ls*(ls+1)*ns/2,1);
    for i=1:ns
        for j=1:ls
            for k=1:j
                %Triangular inferior + diagonal de D, e->d, c->b.
                d((i-1)*ls*(ls+1)/2+j*(j-1)/2+k)=D(j,k+(i-1)*ls);
            end
        end
    end
    
    %Direção nos centros e raio.
    d=[sdpvar(locr,1);d];
    %Restricoes linearizadas.
    SD=S+D;
    %direção viavel simbolica. h(x+d)=h(x)
    linCons=[J*d==zeros(length(J(:,1)),1); d(locr)<=0; r+d(locr)>=0];%; d(1)==0];
    %Restrição linear.
    for j=1:ns
        linCons=linCons+[SD(:,((j-1)*ls+1):(j*ls))>=zeros(ls,ls)];
    end
    
    %-----------------------------------
    
    errc=eps+1;
    errr=eps+1;
    I=0;
    
    %Irrestaurados.
    Vi=replace(V,r,Rp);
    Vi=replace(Vi,c,numc);
    Vi=replace(Vi,S,vS);
    Hi=norm(Vi-b, inf);
    
    %Primeira restauração.
    Rpr=Rp;
    disp(Hi);
    [stat,vSr]=FindS(PosdefCons, S, replace(V-b,c,numc), Hi, r, Rp, 0.1);
    while stat.problem==1
        disp('Primeira restauracao... Hi, Rpr:');
        disp([Hi,Rpr]);
        if Rpr<=eps
            Rpr=Rpr+1;
        else
            Rpr=Rpr*2;
        end
        [stat,vSr]=FindS(PosdefCons, S, replace(V-b,c,numc), Hi, r, Rp, 0.1);
    end
    
    %Matricizando os centros.
    dcmat=sdpvar(lc,cc);
    for j=1:cc
        for k=1:lc
            dcmat(k,j)=d(lc*(j-1)+k);
        end
    end
    
    sdr=d(locr);
    
    while I<=50%errr>eps | errc>eps
        %Direcao tangente.
        [vd, constat]=tanDirYalmip(linCons, d, c, r, S, numc, Rpr, vSr, locr);
        
        if constat.problem==4
            disp('Deu num prob na iteracao:');
            disp(I);
        end
        
        %Se der zero.
        if abs(vd(locr))<eps*eps
            disp('A direção de descida é zero, provavelmente o ponto atual é ótimo.');
            break;
        end
        
        %Se não der pra resolver:
        if constat.problem==1
            disp('SDP linearizado infeasible.');
            [restat,vSr]=iscoveryalmip(PosdefCons+[V==b],S,c,r,numc,Rpr);
            if restat.problem==1
                disp('Não deu para melhorar');
            end
        else
            %Direção no centro.
            vdcmat=replace(dcmat, d, vd);

            %Restaurados.
            Vc=replace(V,c,numc);
            Vr=replace(Vc,r,Rpr);
            Vr=replace(Vr,S,vSr);
            Hr=norm(Vr-b,inf);
            
            %Prints.
            disp('Iteração:');
            disp(I);
            disp('Hi e Hr:');
            disp([Hi,Hr]);
            
            %Calculando erros.
            errrant=errr;
            errcant=errc;
            errc=numc;
            errr=Rpr;
            
            %Direções em S e em r.
            vdsmat=replace(D, d, vd);
            dr=replace(sdr, d, vd);
            
            %Parametros da busca linear.
            buslin=1;
            
            if dalpha==1
                %Tamanhos de passo diferentes.
                for i=1:cc
                    alphac(i)=1/norm(vdcmat(:,i));
                end
                if abs(Rpr+dr)<=eps
                    alphar=0.5;
                else
                    alphar=1;
                end
                alphas=(alphar+alphac*ones(cc,1))/(cc+1);
            else
                %Tamanhos de passo iguais.
                alphaux=min([1,10*errcant,10*errrant]);
                if alphaux<eps
                    alphaux=1;
                end
                alphac=ones(1,cc)*alphaux;
                alphar=1*alphaux;
                alphas=1*alphaux;
            end
            
            %---------------Globalização-----------------
            
            %TESTE.
%             vvv=zeros(1,100);
%             uuu=zeros(1,100);
%                 for kkk=1:100
%                     alphax=(kkk/100);
%                     caplot=numc+alphax*vdcmat;
%                     raplot=Rpr+alphax*dr;
%                     saplot=vSr+alphax*vdsmat;
%                     Vdaux=replace(V,c,caplot);
%                     Vdaux=replace(Vdaux,r,raplot);
%                     Vdaux=replace(Vdaux,S,saplot);
%                     vvv(kkk)=norm(Vdaux-b,2);
%                     uuu(kkk)=raplot;
%                     disp(kkk);
%                 end
%                 vvv=[Hr,vvv];
%                 uuu=[Rpr,uuu];
%                 hold on
%                 plot(1:101,vvv);
%                 text(101,vvv(end),sprintf('%d',I));
%                 %plot(1:101,uuu/norm(uuu));
%                 hold off
            
            while buslin==1
                caux=numc+[alphac;alphac].*vdcmat;
                raux=Rpr+alphar*dr;
                saux=vSr+alphas*vdsmat;
                Vd=replace(V,c,caux);
                Vd=replace(Vd,r,raux);
                Vd=replace(Vd,S,saux);
                Hd=norm(Vd-b,inf);
                if abs(Hr-Hd)<=sigma
                    buslin=0;
                    Rp=raux;
                    numc=caux;
                    vS=saux;
                else
                    alphac=alphac/2;
                    alphar=alphar/2;
                    alphas=alphas/2;
                    if max([alphac,alphar,alphas])<=10^(-12)
                        disp('Comecou a ficar estranho, porque o passo está pequeno demais...');
                        break;
                    end
                end
            end
            
            %--------------Restauração---------------
            Hi=Hd;
            Rprsave=Rpr;
            Rpr=Rp;
            [stat,vSr]=FindS(PosdefCons, S, replace(V-b,c,numc), Hi, r, Rpr, min(0.5,errrant));
            Rpaux=abs(Rpr-Rprsave)/2;
            if Rpaux<=eps*eps;
                Rpaux=0.5;
            end
            restcount=0;
            while stat.problem==1
                if restcount>10
                    Rpr=Rpr+max(0.5, restcount*Rpaux);
                else
                    Rpr=Rpr+Rpaux;
                end
                disp('Em restauracao... Rpr:');
                disp(Rpr);
                [stat,vSr]=FindS(PosdefCons, S, replace(V-b,c,numc), Hi, r, Rpr, min(0.5,errrant));
                restcount=restcount+1;
            end
            
            %Atualiza o erro: |Rpr_{k-1}-Rpr_k|
            errc=norm(numc-errc);
            errr=abs(Rpr-errr);
            
            disp('Raio, centro:');
            disp(Rpr);
            disp(numc);
            disp('Erro centro, erro raio:');
            disp(errc);
            disp(errr);

            %Print dos círculos.
            hold on
            %circle([0;0],1);
            text(0,0,sprintf('V'));
            for j=1:cc
                circlecolor(numc(:,j),sqrt(Rpr),j,I+2);
                text(numc(1,j),numc(2,j),sprintf('%d',I));
            end
            hold off
            
            I=I+1;
            
            disp('----------------------------------');
        end
    end
    disp('Fim.');
    Ropt=Rpr;
    copt=numc;
end