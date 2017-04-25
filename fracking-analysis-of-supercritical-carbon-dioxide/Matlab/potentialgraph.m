function potentialgraph_Jeremy
T1 = [277 283 290 295 300 304];
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
P1_guess = 4*10^6;
for T = [277 290 304.25 315 331]
p = 0.001*10^4:0.001*10^4:2.15*10^4;
u = -2*a*p+(R*T)./(1-b*p)-R*T*log((1-p*b)./p)
plot(p,u,'-r');
xlim([.01*10^4 2.3*10^4])
hold on
end
P = @(V,T) (((R*T)./(V-b)) - (a./(V.^2)));
for i = 1:length(T1)
    P_bin(i) = fsolve(@(P) myfun(P,T1(i)),P1_guess);
    Fs = [P_bin(i) -(P_bin(i)*b+R*T1(i)) a -a*b];
    Vs = (roots(Fs));
    Vs = abs(Vs);
    Vl(i) = min(Vs);
    Vg(i) = max(Vs);
    %plot(Vl(i),P_bin(i),'o',Vg(i),P(Vg(i),T1(i)),'o')
    hold on
end
    
    for i = 1:(length(T1))
    x(i) = Vl(i);
    y(i) = P_bin(i);
    end
    j = 1;
    for i = ((length(T1))+1):(2*length(T1))
    x(i) = Vg((length(T1))-j+1);
    y(i) = P_bin((length(T1))-j+1);
    j = j+1;
    end
    Tnew = [277 283 290 295 300 304];
    rho = 1./x;
    for i = 1:length(Tnew)
    u_bin(i) = -2*a*rho(i)+(R*Tnew(i))./(1-b*rho(i))-R*Tnew(i)*log((1-rho(i)*b)./rho(i));
    
    hold on
    end
    plot([rho(1:6) 7.81*10^3 rho(7:12)],[u_bin 2.18*10^4 flip(u_bin)],'-b');
    i = 1;
for T = [277 283 290 295 300 304];
    
    V = 5*10^-5:0.005*10^-4:5*10^-4;
   % P = @(V,T) (((R*T)./(V-b)) - (a./(V.^2)));
    F_spinodal = [R*T -2*a 4*a*b -2*a*b^2];
    V_ans(:,i) = roots(F_spinodal)
    Vg_sp(i) = V_ans(1,i);
    Vl_sp(i) = V_ans(2,i);
    %plot(V_ans(1,i),P(V_ans(1,i),T),'o',V_ans(2,i),P(V_ans(2,i),T),'o');
    i = i+1;
end
rho_sp_l = 1./Vl_sp
rho_sp_g = 1./Vg_sp
i = 1;
for T = [277 283 290 295 300 304];
ul(i) = -2*a*rho_sp_l(i)+(R*T)./(1-b*rho_sp_l(i))-R*T*log((1-rho_sp_l(i)*b)./rho_sp_l(i));
ug(i) = -2*a*rho_sp_g(i)+(R*T)./(1-b*rho_sp_g(i))-R*T*log((1-rho_sp_g(i)*b)./rho_sp_g(i));   
hold on
i = i+1;
end
plot([rho_sp_g 7.81*10^3 flip(rho_sp_l)],[ug 2.18*10^4 flip(ul)],'-k')

hold off
title('Chemical Potential vs. Density of CO2')
xlabel('Density (p) in mol/m3)');
ylabel('Chemical Potential (J/mol)');
end
function F = myfun(P,T)
a = 0.3649;
b = 4.276*10^-5;
R = 8.314;
Fs = [P -(P*b+R*T) a -a*b];
Vs = (roots(Fs));
Vl = min(Vs);
Vg = max(Vs);
F = R*T*log((Vg-b)/(Vl-b)) + a*((1/Vg) - (1/Vl)) - P*(Vg - Vl);
end
