function [u,v]=flowfield_rankine(R,U,x,y)

% This function provides the (u,v) velocity components for the flow around
% an axisymmetric Rankine body. Note that the x-axis is aligned along the
% axis of the body. 
% Body parameters
% R = Radius of the body
% U = freestream velocity
% x = horizontal position (axis aligned leading to trailing edge)
% y = vertical position

theta = atan2(y,x);
r = sqrt(x.^2+y.^2);
Vr = U*(cos(theta)+R^2./(4*r.^2));
Vt = -U*sin(theta);

u = Vr.*cos(theta) - Vt.*sin(theta);
v = Vr.*sin(theta) + Vt.*cos(theta);

end