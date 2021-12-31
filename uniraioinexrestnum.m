function [Ropt,copt,I,M]=uniraioinexrestnum(PosdefCons, H, S, locr, c, r, numc, Rp, vS, eps, s, symm, mdata, sigma)           
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
    linCons=[-d(locr)>=0];
    if symm==1
        linCons=linCons+[d(1)==0];
    end
    %Restrição linear.
    for j=1:ns
        linCons=linCons+[SD(:,((j-1)*ls+1):(j*ls))>=zeros(ls,ls)];
    end
    
    %-----------------------------------
    
    errc=ones(1,cc)*(eps+1);
    errr=eps+1;
    I=0;
    
    %Irrestaurados.
    Hi=replace(H,r,Rp);
    Hi=replace(Hi,c,numc);
    Hi=replace(Hi,S,vS);
    Hi=norm(Hi, inf);
 
    %Primeira restauração.
    Rpr=Rp;
    disp(Hi);
    [stat,vSr]=FindS(PosdefCons, S, replace(H,c,numc), Hi, r, Rpr, s);
    while stat.problem==1
        disp('Primeira restauracao... Hi, Rpr:');
        disp([Hi,Rpr]);
        if Rpr<=eps
            Rpr=Rpr+1;
        else
            Rpr=Rpr*2;
        end
        [stat,vSr]=FindS(PosdefCons, S, replace(H,c,numc), Hi, r, Rpr, s);
    end
    
    %Matricizando os centros.
    dcmat=sdpvar(lc,cc);
    for j=1:cc
        for k=1:lc
            dcmat(k,j)=d(lc*(j-1)+k);
        end
    end
    
    sdr=d(locr);
    theta=0.999;
    
    datait={0};
    datar=[Rp];
    datarest=[Rpr];
    datac=[(abs(numc.')*ones(length(c(:,1)),1)).'];
    datans=[norm(vS)];
    datansrest=[norm(vSr)];
    datadifs=[0];
    datadr=[0];
    datandc=[zeros(1,cc)];
    datands=[0];
    datalphar=[0];
    datalphac=[zeros(1,cc)];
    datalphas=[0];
    datanh=[Hi];
    datanhrest=[norm(avalcrs(H,c,r,S,numc,Rpr,vSr))];
    datatheta=[theta];
    datas=[s];
    dataerrc=[errc];
    dataerrr=[errr];
    
    %Inserir smax, definir a sequeência de viabilidade previamente e
    %calcular o quanto eu preciso de Hr em função de Hi e s seria a razão entre esses dois.
    
    while  I<=2500% | (errr>eps | max(errc)>eps)
        %Direcao tangente.
        J=jacobnumanual(H,S,c,r,vSr,numc,Rpr,locr,10^(-6)); 
        [vd, constat]=tanDirYalmip(linCons+[J*d==zeros(length(J(:,1)),1)], d, c, r, S, D, dcmat, numc, Rpr, vSr, locr, sigma);
        disp(vd(1:locr));
        %disp('sigma:');
        %disp(sigma*(1-1/(I+2)));
        
        %Decidir direito quem é sigma!
        if constat.problem==0
            disp('Deu success na iteracao:');
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
            [restat,vSr]=iscoveryalmip(PosdefCons+[H==zeros(length(H),1)],S,c,r,numc,Rpr);
            if restat.problem==1
                disp('Não deu para melhorar');
            end
        else
            %Direção no centro.
            vdcmat=replace(dcmat, d, vd);
            disp('vdcmat:');
            disp(vdcmat);
            
            %Restaurados.
            Hr=replace(H,c,numc);
            Hr=replace(Hr,r,Rpr);
            Hr=replace(Hr,S,vSr);
            Hr=norm(Hr,inf);
            
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
           
            %Passo.
            alphaux=1;
            alphac=ones(1,cc)*alphaux;
            alphar=alphaux;
            alphas=alphaux;
            
            %Decidindo s.
            %s=max(min(max(s,1-(norm(errcant)+errrant)/2-10/I), Hr),1);
            s=max(s,1-(norm(errcant)+errrant)/2-1/I);
            %s=min(s,1/(I*I));
            ss=1-s;
            disp('Parametro s:');
            disp(s);
            
            %Decidindo theta.
            theta=min(theta,((s+1)*(Hi-Hr))/(2*(Rpr-Rp-Hr+Hi)));
            
            %Dados.
            if mod(I+1,mdata)==0
                datadr=[datadr;dr];
                datandc=[datandc;(abs(numc.')*ones(length(c(:,1)),1)).'];
                datands=[datands;norm(vdsmat)];
                datatheta=[datatheta;theta];
                datas=[datas;ss];
                datanhrest=[datanhrest;Hr];
            end
            
            %---------------Globalização-----------------
            
            while buslin==1
                caux=numc+[alphac;alphac].*vdcmat;
                raux=Rpr+alphar*dr;
                saux=vSr+alphas*vdsmat;
                Hd=replace(H,c,caux);
                Hd=replace(Hd,r,raux);
                Hd=replace(Hd,S,saux);
                vHd=Hd;
                Hd=norm(Hd,inf);
                %disp(theta*raux+(1-theta)*Hd-(theta*Rp+(1-theta)*Hi)-(Hr-Hi)*(s)/2);
                %hold on
                %plot(vvs,theta*raux+(1-theta)*Hd-(theta*Rp+(1-theta)*Hi)-(Hr-Hi)*(1-vvs)/2);
                %hold off
                if theta*raux+(1-theta)*Hd-(theta*Rp+(1-theta)*Hi)-(Hr-Hi)*(1-s)/2<=0
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
            disp('alpha:');
            disp(alphar);
            
            %--------------Restauração---------------
            
            %DIR NEWTON
            %pn3=[J(locr,:)'*(r-Rpr)];
            %pn2=[J(1:(locr-1),:)'*(c-numc)];
            %pn1=J'*H;
            %pn=pn1+pn2+pn3
            
            Hi=Hd;
            Rprsave=Rpr;
            Rpr=Rp;
            [stat,vSr]=FindS(PosdefCons, S, replace(H,c,numc), Hi, r, Rpr, 1-ss);
            Rpaux=abs(Rpr-Rprsave)/2;
            if Rpaux<=eps*eps;
                Rpaux=0.5;
            end
            restcount=0;
%             dn=-J'*vHd; dn=Rpaux*dn/norm(dn,inf);
%             vdnmat=replace(dcmat, d, dn);
%             disp('Restauracao nos centros:');
%             disp(vdnmat);
%             disp('Restauracao no raio:');
%             disp(dn(locr));
             while stat.problem==1
%                 if restcount>10
%                     Rpr=Rpr+dn(locr)*restcount;
%                     numc=numc+vdnmat;
%                 else
%                     Rpr=Rpr+dn(locr);
%                     numc=numc+vdnmat;
%                 end
%                 disp('Em restauracao... Rpr:');
%                 disp(Rpr);
%                 [stat,vSr]=FindS(PosdefCons, S, replace(H,c,numc), Hi, r, Rpr, s);
%                 restcount=restcount+1;
                if restcount>1
                    Rpr=Rpr+(2^restcount)*Rpaux;
                else
                    Rpr=Rpr+Rpaux;
                end
                
                disp('Em restauracao... Rpr:');
                disp(Rpr);
                [stat,vSr]=FindS(PosdefCons, S, replace(H,c,numc), Hi, r, Rpr, 1-ss);
                restcount=restcount+1;
            end
            
            %Atualiza o erro: |Rpr_{k-1}-Rpr_k|
            errc=(abs(numc-errc).'*ones(length(c(:,1)),1)).';
            errr=abs(Rpr-errr);
%             errr=abs(Rpr-1);
%             errc=norm(numc);
            
            disp('Raio, centro:');
            disp(Rpr);
            disp(numc);
            disp('Erro centro, erro raio:');
            disp(errc);
            disp(errr);
            %3BLUE1BROWN - WELCH LABS

            %Print dos círculos.
            hold on
            circle([0;0],1);
            text(0,0,sprintf('V'));
            for j=1:cc
                circlecolor(numc(:,j),sqrt(Rpr),j,I+2);
                text(numc(1,j),numc(2,j),sprintf('%d',I));
            end
            hold off
            
            I=I+1;
            
            %Dados 2.
            if mod(I,mdata)==0
                datait=[datait;{I}];
                datar=[datar;Rp];
                datarest=[datarest;Rpr];
                datac=[datac;(abs(numc.')*ones(length(c(:,1)),1)).'];
                datans=[datans;norm(vS)];
                datansrest=[datansrest;norm(vSr)];
                datadifs=[datadifs;norm(vS-vSr,inf)];
                datalphar=[datalphar;alphar];
                datalphac=[datalphac;alphac];
                datalphas=[datalphas;alphas];
                datanh=[datanh;Hi];
                dataerrc=[dataerrc;errc];
                dataerrr=[dataerrr;errr];
            end
            
            disp('----------------------------------');
        end
    end
    disp('Fim.');
    Ropt=Rpr;
    copt=numc;
    
    VarName={'It','r','r_Rest','c','nS','nS_Rest','nH','nH_Rest','dr','dc','ndS','alphar','alphac','alphas','erro_c','erro_r'};
    M=table(datait,datar,datarest,datac,datans,datansrest,datanh,datanhrest,datadr,datandc,datands,datalphar,datalphac,datalphas,dataerrc,dataerrr,'VariableNames',VarName);
end
