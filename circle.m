%Recebe: um centro x e um raio r.
%Retorna: um plot do circulo.

function h = circle(x,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x(1);
yunit = r * sin(th) + x(2);
h = plot(xunit, yunit, 'LineWidth', 1,'Color', [0,0,0]);
%plot(x(1),x(2),'.');
%h = plot(xunit, yunit,'--','Color',[0.7,0.7,0.7]);
hold off
end