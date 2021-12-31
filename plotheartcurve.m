function [X,Y,Z]=plotheartcurve(a,b,h)
    x=a:h:b;
    [X,Y]=meshgrid(x);
    Z=(X.^2+Y.^2-1).^3-(X.^2).*(Y.^3);
    hold on
    contour(X,Y,Z,[0,0],'Color',[0,0,0]);
    hold off
end