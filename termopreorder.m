%Recebe: p=vetor de geradores do preorder. s=coef SOS do termo. str=string
%de potencias dos fatores na ordem: p1,p2,...,pe,q1,...,qm.
%Retorna: o termo T correspondente ao coef s e o expoente exp.

function T=termopreorder(p,s,exp)
    T=s;
    for i=1:length(exp)
        T=T*(p(i)^str2double(exp(i)));
    end
end