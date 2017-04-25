function xmin1=streamline_rankine(R,U,ID)

%This function plots the streamlines and shape of a Rankine body

%domain for plotting streamlines
nx = 800;  xmin = -10*R; xmax = 3.5*R;
ny = 200;  ymin = -2*R;  ymax = 2*R;
[x,y] = meshgrid(linspace(xmin, xmax, nx), ...
linspace (ymin, ymax, ny));

%FLOW OVER CYLINDER
theta = atan2(y,x);
r = sqrt(x.^2 + y.^2);
psi = U/2*( y.^2 - R.^2/2*( 1 + x./sqrt(x.^2+y.^2) )) ;
psimax = max(psi(:));    
levels=[0.15:8]/8*psimax;

%STREAMLINES OF BODY
clf
%contourf(x, y, U*sqrt(1+(R./r).^2/2.*cos(theta)+(R./r).^4/16),[0:50]/50*U*1.3);
hold on
[C,h]=contour(x, y, psi, levels);
set(h,'color',[0.4 0.7 1]);
axis equal
axis([xmin xmax ymin ymax])
xlabel('x')
ylabel('y')
hold on

thetamin=atan2(R,xmax)*0.95;
thetamax=acos(2*ID^2-1);
tt = ([thetamin:0.01:thetamax]);
ytube=R*sqrt((1+cos(tt))/2);
rtube=ytube./sin(tt);

plot([xmin xmax],[0 0],'color',[0.5 0.5 0.5]);
fill([rtube.*cos(tt) xmax],[ytube ID*R],0.7*[1 1 1]);
fill([rtube.*cos(tt) xmax],-[ytube ID*R],0.7*[1 1 1]);
plot([rtube.*cos(tt) xmax],[ytube ID*R],'k');
plot([rtube.*cos(tt) xmax],-[ytube ID*R],'k');
plot(rtube(end)*cos(tt(end))*[1 1], ytube(end)*[-1 1],'k');
%gets xmin
ytube1=R*sqrt((1+cos(thetamax))/2);
rtube1=ytube1./sin(thetamax);
xmin1=rtube1*cos(thetamax);

end