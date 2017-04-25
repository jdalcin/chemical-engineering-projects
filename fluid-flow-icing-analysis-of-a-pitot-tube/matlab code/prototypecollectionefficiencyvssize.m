function [D,CollectionEfficiency]=prototypecollectionefficiencyvssize
%Prototype constants
VFo=250; %freestream velocity of air (m/s)
rhop= 999.97; %water droplet density (kg/m^3)
rhof=.4; %air density (kg/m^3)
u=1.4*10^-5; %air dynamic viscosity (kg/(m*s))
g=-9.8; %acceleration of gravity (m/s^2)
R=5*10^-3; %pitot outer tube radius (m)
Ri=2*10^-3; %pitot inner tube radius (m)
D= linspace(2*10^-6,80*10^-6,100); %potential diameters of water droplet (m)
%Laboratory Constants
% VFo=50; %freestream velocity of air (m/s)
% rhop= 999.97*3; %water droplet density (kg/m^3)
% rhof=1.2; %air density (kg/m^3)
% u=1.711*10^-5; %air dynamic viscosity (kg/(m*s))
% g=-9.8; %acceleration of gravity (m/s^2)
% D=((.4*250*Dp)/(1.4*10^-5))*((u)/(1.2*50))
% Ri=((4*10^-3/Dp)*D)/2 %pitot inner tube radius (m)
% R=Ri./.4 %pitot outer tube radius (m)
% Dp=D;

xmin1=streamline_rankine(R,VFo,Ri/R); %x minimum of pitot tube
options = odeset('InitialStep',1e-4,'MaxStep',1e-4,'RelTol',1e-6,'AbsTol',1e-6);%,'Events',@events); %step from Kiger's Matlab HW2 Code
%comes up with particle streamlines up until ones which hit the inner
%radius of pitot tube for various particle sizes
for k=1:length(D)
    Dp=D(k);
    i=0;
    Yin=Ri;
    mp=rhop*((pi*Dp.^3)/6); %potential masses of water droplet
    while (Yin>=Ri) 
    [VFox,VFoy]=flowfield_rankine(R,VFo,-.1,(Ri*(100-i))/100); %initial air velocities: x and y
    VPox=VFox; %%initial particle velocity wrt air:x-direction
    VPoy=VFoy; %initial particle velocity wrt air: y-direction
    initial=[-.1 VPox (Ri*(100-i))/100 VPoy];
    [T,Coupled]=ode15s(@(t,init)coupled_ode(t,init,Dp,g,rhop,rhof,u,R,VFo,mp),[0 .1/VFo],initial,options);
    [Xin,Yin]=intersections(Coupled(:,1),Coupled(:,3), ones(length(Coupled(:,1)),1)*xmin1,linspace(-2*R,2*R,length(Coupled(:,1)))');
     Yin=Yin(1);
     Y(i+1,1)=Yin;
    i=i+1;
    end
Y; %Y vector of pitot tube intersection points at various streamlines
CollectionEfficiency(k,1)=((Ri*(100-i)/100)^2/(Ri)^2)*100;
end
CollectionEfficiency;
%plots Collection Efficiences vs. Particle Diameter Size
figure(2)
plot(D,CollectionEfficiency)
ylabel('Collection Efficiency')
xlabel('Particle Diameter (m)')
title('Pitot Tube Collection Efficiencies')
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