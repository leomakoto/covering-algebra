function plottrianglecurve(h)
    a=sqrt(3/4);
    x1=-0.5:h:0;
    y1=a+2*a*x1;
    x2=0:h:0.5;
    y2=a-2*a*x2;
    x3=-0.5:h:0.5;
    y3=0.*x3;
    hold on
        plot(x1,y1,'Color',[0,0,0]);
        plot(x2,y2,'Color',[0,0,0]);
        plot(x3,y3,'Color',[0,0,0]);
    hold off
end