%%% Orbital motion simulation: orbital evolution part
%%% - orbital evolution of a rotating asteroid illimuniated by the Sun
%%% Author: Selver Pepic
%%% Last major edit: 17.06.2020.

%%% Main steps:
%%% - takes theoretical values for net thermal recoil forces (because the thermal code did not fully work)
%%% - uses the 2nd order Verlet energy-conserving time stepping algorithm to calculate orbital evolution over few centuries


%% theoretical calculation of Yarkowsky effect
% physical parameters
k = 0.3; % thermal conductivity, 2.5 W/mK or 0.3  
rho = 1190; % density, 2600 kg/m3 or 2000 
cv = 500; % heat capacity, 2000 J/kg?C or 600
sig = 5.67E-8; % stefan boltzmann constant, 5.67E-8
c = 3*10^8; % speed of light 3e8 m/s

Rasteroid = 250; % asteroid radius (Bennu = 493 m, a_distance = 1.10 AU)
% da_dt = ?18.95 � 0.10  (Bennu)
Tday = 4.288 * 60 * 60;
w = 2*pi/Tday;
%w = 0.01; % rad/s, angular velocity, positive = counterclockwise
I = 1000; % W/m^2, solar flux on asteroid (1367 W/m2 at Earth)

% normalized variables (theoretical)
L_diff = sqrt(k/(rho*cv*w));
R_norm = Rasteroid/L_diff;
G = sqrt(k*rho*cv);
T_sub = (I/sig)^(1/4);
T_norm = G*sqrt(w)/(sig*T_sub^3); % (k*T_sub/L_diff) / I;
% Force factors
W = -0.5*T_norm/(1+T_norm+0.5*T_norm^2);
FI = pi*Rasteroid^2 * I / (c * rho*4/3*pi*Rasteroid^3); % I*S/mc = acceleration
gamma = 0; % spin axis angle from the orbital plane normal;
Tyear = 1.20 * 365*24*60*60;
wyear = 2*pi/Tyear;
dadt = -8/9 * FI/wyear * W * cos(gamma); % (m/s^2)/(1/s)) = m/s
% mreasured to be -18.95 * 10^-4;
dadt_auyr = dadt / (1.496*10^11) * 365*24*60*60;
dadt_auMyr = dadt_auyr * 10^6;
a_yark = 4/9 * FI * W;


%% setup
% numerical
tmax = 1*10^2;
niter = 10*10^5;
dt = tmax/niter;
t = dt:dt:tmax;
% physics
g = 1; % gsol at 1 AU = 5.9E-3
gsol = 5.9*10^-3;
Ro2 = 1; % GM = gRo2
a_yark_t = 10^8*a_yark/gsol;%-0.005;
xo = 1.36;
yo = 0;
vox = 0;
voy = sqrt(g*Ro2/xo); % 0.765
Tperiod = xo*2*pi/voy;
% v^2/r = GM/r^2, GM = gRo^2
% v^2 = gRo^2/r
% v^2 = GM = gRo2  stable orbit
% 0.5*v^2 = GM/r = gRo2/r
% v = sqrt(2) * gRo2/r

% iterate
x = zeros(niter,1);
y = zeros(niter,1);
vx = zeros(niter,1);
vy = zeros(niter,1);
x(1) = xo;
y(1) = yo;
vx(1) = vox;
vy(1) = voy;
ax = zeros(niter,1);
ay = zeros(niter,1);
theta = zeros(niter,1);

%% main iteration
for it = 1:niter
    
    theta(it) = atan2(y(it),x(it));
    a = g*Ro2/(x(it)^2+y(it)^2);
    ax(it) = -a*cos(theta(it)) + a_yark_t * (-sin(theta(it)));
    ay(it) = -a*sin(theta(it)) + a_yark_t * cos(theta(it));
    
    % Euler forward
    %vx(it+1) = vx(it) + ax(it)*dt;
    %vy(it+1) = vy(it) + ay(it)*dt;
    %x(it+1) = x(it) + vx(it)*dt;
    %y(it+1) = y(it) + vy(it)*dt;
    
    % symplectic Euler-Cromer (1st order)
    %vx(it+1) = vx(it) + ax(it)*dt;
    %vy(it+1) = vy(it) + ay(it)*dt;
    %x(it+1) = x(it) + vx(it+1)*dt;
    %y(it+1) = y(it) + vy(it+1)*dt;
    
    % verlet
    if it == 1
        %vx(it+1) = vx(it) + ax(it)*dt;
        %vy(it+1) = vy(it) + ay(it)*dt;
        %x(it+1) = x(it) + vx(it+1)*dt;
        %y(it+1) = y(it) + vy(it+1)*dt;        
        x_minus = x(it) - dt*vx(it) + 0.5*ax(it)*dt^2;
        y_minus = y(it) - dt*vy(it) + 0.5*ay(it)*dt^2;
        x(it+1) = 2*x(it) - x_minus + ax(it)*dt^2;
        y(it+1) = 2*y(it) - y_minus + ay(it)*dt^2;
        vx(it) = (x(it+1) - x_minus) / 2/dt;
        vy(it) = (y(it+1) - y_minus) / 2/dt;
    else
        x(it+1) = 2*x(it) - x(it-1) + ax(it)*dt^2;
        y(it+1) = 2*y(it) - y(it-1) + ay(it)*dt^2;
        vx(it) = (x(it+1) - x(it-1)) / 2/dt;
        vy(it) = (y(it+1) - y(it-1)) / 2/dt;
    end
    
    % verlet velocity
    %if it == 1
    %    vx(it+1) = vx(it) + ax(it)*dt;
    %    vy(it+1) = vy(it) + ay(it)*dt;
    %    x(it+1) = x(it) + vx(it+1)*dt;
    %    y(it+1) = y(it) + vy(it+1)*dt;
    %else
    %    x(it+1) = x(it) + vx(it)*dt + 1/2*ax(it)*dt^2;
    %    y(it+1) = y(it) + vy(it)*dt + 1/2*ay(it)*dt^2;
    %        theta(it+1) = atan2(y(it+1),x(it+1));
    %        a = g*Ro2/(x(it+1)^2+y(it+1)^2);
    %        ax(it+1) = -a*cos(theta(it+1));
    %        ay(it+1) = -a*sin(theta(it+1));
    %    vx(it+1) = vx(it) + 0.5*dt*(ax(it)+ax(it+1));
    %    vy(it+1) = vy(it) + 0.5*dt*(ay(it)+ay(it+1));
    %end
    
end
x = x(1:niter);
y = y(1:niter);
vx = vx(1:niter);
vy = vy(1:niter);
theta = theta(1:niter);

%% figures
figure(1)
    % Earth
    ang = linspace(0,2*pi,100);
    plot(cos(ang),sin(ang),'r-');
    hold on;
    scatter(0,0,'y','O','LineWidth',5);
    % asteroid
    hold on;
    scatter(x,y,'b.');
    hold on;
    scatter(xo,yo,'b','O','LineWidth',2);
    hold on;
    scatter(x(end),y(end),'b','X','LineWidth',2);
    xlabel('x (a.u.)')
    ylabel('y (a.u.)')
    title('Asteroid�s orbit vs. time')
    %xlim([-1,1])
    %ylim([-1,1])
%figure(2)
%    plot(t,x,'b');
%    hold on;
%    plot(t,y,'r')  
%    title('x and y vs t')
figure(3)
    plot(t,sqrt(x.^2+y.^2))
    title('distance vs t')
figure(4)
    t_period = 0;%zeros(round(niter/(Tperiod/dt))-2,1);
    x_period = 0;%zeros(round(niter/(Tperiod/dt))-2,1);
    y_period = 0;%zeros(round(niter/(Tperiod/dt))-2,1);
    count = 0;
    for it = 1:round(Tperiod/dt):niter
        if it*round((Tperiod/dt)) >= niter
            break;
        end
        count = count + 1;
        t_period(count) = t(it);
        x_period(count) = x(it);
        y_period(count) = y(it);        
    end
    plot(t_period,sqrt(x_period.^2+y_period.^2))
    title('distance vs t')
figure(5)
    E = 1/2*(vx.^2+vy.^2) - g*Ro2./sqrt(x.^2+y.^2);
    plot(t,E)
    E_error = (max(E)-min(E)) / mean(E);
    title(['E vs t, error = ',num2str(E_error)])
    
    %% test
    t1 = t(1:10*round(Tperiod/dt));
    x1 = x(1:10*round(Tperiod/dt));
    y1 = y(1:10*round(Tperiod/dt));
    d1 = sqrt(x1.^2+y1.^2);
    figure(1)
        scatter(x1,y1,'b.');
    figure(2)
        plot(t1,d1)
        title('distance vs t')
        
    t2 = t(niter-10*round(Tperiod/dt):niter);
    x2 = x(niter-10*round(Tperiod/dt):niter);
    y2 = y(niter-10*round(Tperiod/dt):niter);
    d2 = sqrt(x2.^2+y2.^2);
    figure(3)
        scatter(x2,y2,'b.');
    figure(4)
        plot(t2,d2)
        title('distance vs t')
        
      %% test 2
      r = sqrt(x.^2 + y.^2);
      k_mean = 100*round(Tperiod/dt);
      r_movmean = movmean(r,k_mean); 
      figure(5)
        plot(t,r_movmean-r_movmean(1))
        title('Semi-major distance drift vs. time')
        xlabel('t (a.u.')
        ylabel('\Delta r (a.u.)')
      