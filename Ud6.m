function [loc] = Ud6( radius )
x1=3*radius*rand();
y1=radius*sqrt(3)/2*rand();
loc=zeros(1,2);
if  y1<=(sqrt(3)*radius-x1*sqrt(3))
    x=x1;
    y=y1;
elseif  y1<=(-sqrt(3)*2*radius+x1*sqrt(3))
    x=x1-3*radius;
    y=y1;
else
    x=x1-1.5*radius;
    y=y1-sqrt(3)*radius/2;
end
loc=[x,y] ;
end

