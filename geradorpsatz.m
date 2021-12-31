%%---------------Dados de entrada-----------------.

n=2; %Dimensao.
e=3; %Numero de bolas.
d=2; %Metade do grau dos SOS.
b=1; %Metade do grau do termo do monoide, para Stengle.
psatzver='putinar'; %Versao do Psatz.
var=sdpvar(n,1); %Variaveis dos polinomios.
dq=2; %Vetor coluna com os graus dos polinomios de V.
q=[-var(1)^2-var(2)^2+1.0]; %Polinomios que definem V com \geq em coluna.
m=length(q); %Quantidade de polinomios que definem V.

%%-----------------Coisas comuns aos dois Psatze-----------------.

if strcmp(psatzver,'putinar')
    ns=1+e+m; %Quantidade de variaveis SOS de Putinar.
elseif strcmp(psatzver,'stengle')
    ns=2^(e+m); %Quantidade de variaveis SOS de Stengle.
else
    disp('Ainda não utilizo essa versão do P-satz.');
end

ls=nchoosek(n+d,d); %Tamanho de cada variavel SOS.
locr=n*e+1; %Posicao do raio.
x=sdpvar(locr+ls*(ls+1)*ns/2,1); %Vetorzao de variaveis (renaming).
c=sdpvar(n,e,'full');

%Criando os centros.
for j=1:e
    for k=1:n
        c(k,j)=x(n*(j-1)+k);
    end
end

%Criando o raio.
r=x(locr);

%Criacao dos circulos.
p=sdpvar(e,1);
for i=1:e
    p(i)=(var-c(:,i)).'*(var-c(:,i))-r;
end

v=monvecyalmip(var,d); %Vetor de monomios [x]_d.

%Definindo os polinomios SOS e o polinomio Psatz.
S=sdpvar(0,0);
for i=1:ns
    S=[S,sdpvar(ls,ls,'symmetric','real')];
end

%Vetorizando os coeffs SOS.
for i=1:ns
    for j=1:ls
        for k=1:j
            %Triangular inferior + diagonal de D, e->d, c->b.
            x(locr+((i-1)*ls*(ls+1)/2+j*(j-1)/2+k))=S(j,k+(i-1)*ls);
        end
    end
end

%%-----------------Geracao do Psatz-----------------.

if strcmp(psatzver,'putinar')
    disp('Psatz de Putinar');
    
    polvet=[1;p;q]; %Vetor total de polinomios do problema.

    %Calculando Psatz de Putinar.
    Psatz=1;
    for i=1:ns
        Psatz=Psatz+polvet(i)*(v.'*S(:,((i-1)*ls+1):(i*ls))*v);
    end
    
    %Extraindo os coeficientes.
    H=coefficients(Psatz, var);

elseif strcmp(psatzver,'stengle')
    disp('Psatz de Stengle');
    
    exp=dec2bin(0:(ns-1));
    polvet=[p;q];
    
    %Calculando Psatz de Stengle.
    Psatz=0;
    for i=1:ns
        Psatz=Psatz+termopreorder(polvet,(v.'*S(:,((i-1)*ls+1):(i*ls))*v),exp(i,:));
    end
    
    %Termo do monoide mult.
    mono=1;
    for i=1:e
       mono=mono*p(i)^(2*b);
    end
    Psatz=Psatz+mono;
    
    %Extraindo os coeficientes.
    H=coefficients(Psatz, var);
else
    disp('Ainda não utilizo essa versão do P-satz.');
end

%%---------------Texto---------------

psatxt=sdisplay(H); %Em forma de symcell - transformar cada uma para char.
nh=length(psatxt);
fname=sprintf('%s_n=%d_e=%d_d=%d.txt',psatzver,n,e,d);
fid=fopen(fname,'w');
fprintf(fid,'function h(x)\r\n     real, dimension(%d) :: x\r\n     real, dimension(%d) :: h\r\n\r\n', locr+ls*(ls+1)*ns/2, nh);
for i=1:nh
    auxtext=char(psatxt(i));
    auxtext=strrep(auxtext,'^','**');
    fprintf(fid,'     h(%d) = %s\r\n',i,auxtext);
end
fprintf(fid,'end function');
fclose(fid);

%%---------------Jacobiana do Psatz para Putinar---------------

M=v*v';
tr=sdpvar(ns,1);
for i=1:ns
    tr(i)=trace(M*S(:,((i-1)*ls+1):(i*ls)));
end
J=sdpvar(1,locr);

%Derivadas nos centros.
for i=1:e
    for j=1:n
        J(n*(i-1)+j)=tr(i+1)*2*(c(j,i)-var(j));
    end
end

%Derivada no raio.
J(locr)=-[0,ones(1,ns-m-1),zeros(1,m)]*tr;
Js=sdpvar(ls*(ls+1)*ns/2,1).';

%Derivadas nos S.
for i=1:ns
    for j=1:ls
        for k=1:j
            %Triangular inferior + diagonal, e->d, c->b.
            aux=polvet(i)*M(j,k);
            if j~=k
                aux=aux*2;
            end
            Js((i-1)*ls*(ls+1)/2+j*(j-1)/2+k)=aux;
        end
    end
end

J=[J,Js]; %Jacobiano total em forma de polinomio.

%Extracao de coeficientes de cada componente de J.
dp=2*d+max([2;dq]); %Grau do Psatz.
Jtotal=sdpvar(nchoosek(n+dp,dp),length(J));

%Reordenacao de mvec na mesma ordem que "coefficients" retorna.
mvec=monvecyalmip(var,dp);
nmvec=length(mvec);
[~,mvec]=coefficients(ones(1,nmvec)*mvec, var);

%Montagem do Jacobiano total.
for i=1:length(J)
    [auxc,auxv]=coefficients(J(i),var);
    for j=1:nmvec
        %Para depois q acabar o auxv.
        if length(auxv)<j
            auxv=vertcat(auxv,mvec(j));
            auxc=vertcat(auxc,0);
        %Para encaixar o zero onde nao tem, pois "coefficients" retorna uma coisa sem zeros.
        elseif ~isequal(auxv(j),mvec(j))
            auxv=vertcat(auxv(1:(j-1)),mvec(j),auxv(j:end));
            auxc=vertcat(auxc(1:(j-1)),0,auxc(j:end));
        end
    end
    Jtotal(:,i)=auxc;
end

%%---------------Texto---------------

jacobtxt=sdisplay(Jtotal); %Em forma de symcell - transformar cada uma para char.
njac=length(jacobtxt);
fname=sprintf('putinar_n=%d_e=%d_d=%d_jacobiana.txt',n,e,d);
fid=fopen(fname,'w');
ljac=length(Jtotal(:,1));
cjac=length(Jtotal(1,:));
fprintf(fid,'function jacobh(x)\r\n     real, dimension(%d) :: x\r\n     real, dimension(%d,%d) :: jacobh\r\n\r\n', locr+ls*(ls+1)*ns/2, ljac, cjac);
for i=1:ljac
    for j=1:cjac
        auxtext=char(jacobtxt(i,j));
        auxtext=strrep(auxtext,'^','**');
        fprintf(fid,'     jacobh(%d,%d) = %s\r\n',i,j,auxtext);
    end
end
fprintf(fid,'end function');
fclose(fid);