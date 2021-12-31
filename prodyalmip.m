%Produtorio simbolico. Recebe um vetor e retorna o produto de seus
%elementos simbolicos.
%Obs: complemento a monvecplus.

function p=prodyalmip(v)
    p=1;
    for i=1:length(v)
        p=p*v(i);
    end
end