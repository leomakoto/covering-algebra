%Recebe: um centro x e um raio r.
%Retorna: um plot do circulo.

function h = circlecolor(x,r,i,j)
if j==0
    col=zeros(1,3);
else
    col=ones(1,3)*0.7;
     i=mod(i,3)+1;
     col(i)=1/j;
    %Melhorar.
end
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x(1);
yunit = r * sin(th) + x(2);
h = plot(xunit, yunit ,'Color', col);
%plot(x(1),x(2),'.');
hold off
end