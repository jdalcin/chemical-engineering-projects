function collectionefficiencygravity
%Prototype constants
VFo=250; %freestream velocity of air (m/s)
rhop= 999.97; %water droplet density (kg/m^3)
rhof=.4; %air density (kg/m^3)
u=1.4*10^-5; %air dynamic viscosity (kg/(m*s))
g=-9.8; %acceleration of gravity (m/s^2)
R=5*10^-3; %pitot outer tube radius (m)
Ri=2*10^-3; %pitot inner tube radius (m)
D= linspace(2*10^-6,80*10^-6,100); %potential diameters of water droplet (m)
 Dp=80*10^-6;
%Laboratory Constants
% VFo=50; %freestream velocity of air (m/s)
% rhop= 999.97*3; %water droplet density (kg/m^3)
% rhof=1.2; %air density (kg/m^3)
% u=1.711*10^-5; %air dynamic viscosity (kg/(m*s))
% g=-9.8; %acceleration of gravity (m/s^2)
% D=((.4*250*Dp)/(1.4*10^-5))*((u)/(1.2*50));
% Ri=((4*10^-3/Dp)*D)/2; %pitot inner tube radius (m)
% R=Ri./.4; %pitot outer tube radius (m)
% Dp=D;
 mp=rhop*((pi*Dp.^3)/6); %potential masses of water droplet
%initial particle velocities
i=0;
Yin=Ri;
xmin1=streamline_rankine(R,VFo,Ri/R); %x minimum of pitot tube
options = odeset('InitialStep',1e-4,'MaxStep',1e-4,'RelTol',1e-6,'AbsTol',1e-6);; %step from Kiger's Matlab HW2 Code
%comes up with particle streamlines up until ones which hit the inner
%radius of pitot tube
while (Yin>=Ri) 
    [VFox,VFoy]=flowfield_rankine(R,VFo,-.1,(Ri*(100-i))/100); %initial air velocities: x and y
    VPox=VFox; %%initial particle velocity wrt air:x-direction
    VPoy=(mp*g)/(3*pi*u*Dp)+VFoy; %initial particle velocity wrt air: y-direction
    initial=[-.1 VPox (Ri*(100-i))/100 VPoy];
    [T,Coupled]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .1/VFo],initial,options);
    [Xin,Yin]=intersections(Coupled(:,1),Coupled(:,3), ones(length(Coupled(:,1)),1)*xmin1,linspace(-2*R,2*R,length(Coupled(:,1)))');
     Yin=Yin(1);
     Y(i+1,1)=Yin;
    i=i+1;
end
initial=[-.1 VPox -(Ri*(100-i))/100 VPoy];
     [T2,Coupled2]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .1/VFo],initial,options);
     [uf2,vf2]=flowfield_rankine(R,VFo,Coupled2(:,1),Coupled2(:,3));
initial=[-.1 VPox -(1.25*R) VPoy];
     [T3,Coupled3]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .01],initial,options);
initial=[-.1 VPox (1.25*R) VPoy];
     [T4,Coupled4]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .01],initial,options);
initial=[-.1 VPox -(R) VPoy];
     [T3,Coupled5]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .01],initial,options);
initial=[-.1 VPox (R) VPoy];
     [T4,Coupled6]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .01],initial,options);
%-----
 Y; %Y vector of pitot tube intersection points at various streamlines
 CollectionEfficiency=((Ri*(100-(i-1))/100)^2/(Ri)^2)*100
% initial=[-70*R VPox .005 VPoy];
%     [T3,Coupled3]=ode45(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,.005,VFo,mp),[0 .01],initial);
% initial=[-70*R VPox 1.5179*10^-3 VPoy];
%     [T4,Coupled]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,.005,VFo,mp),[0 .01],initial);

%calculates Drag Force, Reynolds Number, particle relative velocity, and
%fluid and particle velocity using general drag

% [uf,vf]=flowfield_rankine(R,VFo,Coupled(:,1),Coupled(:,3));
% Vx=Coupled(:,2);
% Vy=Coupled(:,4);
% sqrtr4=sqrt((Vx-uf).^2+(Vy-vf).^2);
% Reb4=(rhof*Dp*sqrtr4)./u; %Reynold's number
% Cd3=@(Rep)(24./Rep).*(1+.27*Rep).^.43+.47*(1-exp(-.04*Rep.^.38));
% Cdd=Cd3(Reb4);
% k=@(Cd)(rhof*Cd*pi*Dp^2)/(8*mp);
% FDdx=@(Vx,Vy,Vfx,Vfy)mp*k(Cdd).*(uf-Vx).*sqrt((uf-Vx).^2+(vf-Vy).^2);
% b=FDdx(Vx,Vy,uf,vf);
% FDdy=@(Vx,Vy,Vfx,Vfy)mp*k(Cdd).*(vf-Vy).*sqrt((uf-Vx).^2+(vf-Vy).^2);
% FDd=sqrt(FDdx(Vx,Vy,uf,vf).^2+FDdy(Vx,Vy,uf,vf).^2);

% calculates Drag Force, Reynolds Number, particle relative velocity, and
%fluid and particle velocity using modified stokes drag

[uf,vf]=flowfield_rankine(R,VFo,Coupled(:,1),Coupled(:,3));
Vx=Coupled(:,2);
Vy=Coupled(:,4);
sqrtr4=sqrt((Vx-uf).^2+(Vy-vf).^2);
Reb4=@()(rhof*Dp*sqrtr4)./u;
z=(1+Reb4()./(4*(1+sqrt(Reb4())))+Reb4()./60);
FDdx=@(Vx,uf)mp*((18*u)/(rhop*Dp^2))*(uf-Vx).*z;
FDdy=@(Vy,vf)mp*((18*u)/(rhop*Dp^2))*(vf-Vy).*z;
FDd=sqrt(FDdx(Vx,uf).^2+FDdy(Vy,vf).^2);


%plots streamlines
figure(1)
hold all
streamline_rankine(R,VFo,Ri/R);
plot(Coupled(:,1),Coupled(:,3),'-r')
plot(Coupled2(:,1),Coupled2(:,3),'-r')
plot(Coupled3(:,1),Coupled3(:,3),'-r')
plot(Coupled4(:,1),Coupled4(:,3),'-r')
plot(Coupled5(:,1),Coupled5(:,3),'-r')
plot(Coupled6(:,1),Coupled6(:,3),'-r')
% plot(ones(length(Coupled(:,1)),1)*xmin1,linspace(0,2*R,length(Coupled(:,1)))','-g')
title('Trajectory of Particle in Air')
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
ylabel('Vertical Distance From Center of Pitot Tube (m)')
hold off all
%plots Reynolds number
figure('Name','Reynolds','NumberTitle','Off')
semilogy(Coupled(:,1),Reb4())
 xlim([-.1 xmin1]);
ylim([10^-.5 150])
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
ylabel('Reynolds Number')
title('Particle Reynolds Number over Various Distances from Pitot Tube')
%plots drag force & particle weight
figure(4)
s=length(Coupled(:,1));
wp=-mp*g*ones(s,1);
semilogy(Coupled(:,1),FDd,'-b',Coupled(:,1),wp,'-r')
xlim([-.1 xmin1]);
title('Drag Force Variation')
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
ylabel('Drag Force (N)')
legend('Drag Force','Particle Weight','Location','BestOutside')
%plots particle relative velocity
figure(5)
semilogy(Coupled(:,1),sqrtr4)
xlim([-.1 xmin1]);
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
ylabel('Particle Relative Velocity (m/s)')
title('Particle Relative Velocity Variation')
% %plots fluid and particle velocity
hold all
figure(6)
fp=[uf vf Vx Vy];
plot(Coupled(:,1),fp)
title('Fluid & Particle Velocities')
ylabel('Velocity (m/s)')
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
legend('Uf','Vf','Up','Vp','Location','BestOutside')
 ylim([-10 260])
 xlim([-.1 xmin1-.0005]);
hold all off
figure(7)
hold all
fp=[uf2 vf2 Coupled2(:,2) Coupled2(:,4)];
plot(Coupled2(:,1),fp)
title('Fluid & Particle Velocities')
ylabel('Velocity (m/s)')
xlabel('Horizontal Distance From Center of Pitot Tube (m)')
legend('Uf','Vf','Up','Vp','Location','BestOutside')
 ylim([-200 260])
 xlim([-.1 xmin1-.0005]);
end

function D=coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp)
% air velocities
[VFx,VFy]=flowfield_rankine(R,VFo,init(1),init(3));
%Reynold's Number of particle
sqrtre=sqrt((init(2)-VFx).^2+(init(4)-VFy).^2);
k=@(Cd)(rhof*Cd*pi*Dp^2)/(8*mp);
Rep=@()(rhof*Dp*sqrtre)./u;
%odes using stokes drag
% D(1)=init(2); %ode of particle x-velocity
% D(2)=((18*u)/(rhop*Dp^2))*(VFx-init(2)); %ode of particle x-acceleration
% D(3)=init(4); %ode of particle y-velocty
% D(4)=((18*u)/(rhop*Dp^2))*(VFy-init(4))+g; %ode of particle y-acceleration
%odes using general drag
% Cd=@(Rep)(24./Rep).*(1+.27*Rep).^.43+.47*(1-exp(-.04*Rep.^.38));
% D(1)=init(2);
% D(2)=k(Cd(Rep())).*(VFx-init(2)).*sqrt((VFx-init(2)).^2+(VFy-init(4)).^2);
% D(3)=init(4);
% D(4)=k(Cd(Rep())).*(VFy-init(4)).*sqrt((VFx-init(2)).^2+(VFy-init(4)).^2)+g;
%odes using modified stokes drag
z=(1+Rep()/(4*(1+sqrt(Rep())))+Rep()/60);
D(1)=init(2); %ode of particle x-velocity
D(2)=((18*u)/(rhop*Dp^2))*(VFx-init(2))*z; %ode of particle x-acceleration
D(3)=init(4); %ode of particle y-velocty
D(4)=((18*u)/(rhop*Dp^2))*(VFy-init(4))*z+g; %ode of particle y-acceleration
D=D';
end