%%% Yarkowsky effect simulation: thermal/radiative part
%%%  - temperature map of a rotating asteroid illimuniated by the Sun, and the net thermal recoil forces
%%% Author: Selver Pepic
%%% Last major edit: 17.06.2020.

%%% Main steps:
%%% - radiative transfer from the "Sun" via temperature gradient boundary conditions
%%% - FDM heat diffusion to calculate the temperature field
%%% - calculation of black body radiation and resulting thermal photon recoil
%%% - calcualation of net thermal photon recoil, decomposition into tangential (along the orbit) and normal parts

%%% Note:
%%% numerical scheme has an unclear error that makes it work only for illumination from 0 and 90 degrees
%%% in other cases, there is a clear disbalance in energy input and thermal energy increase
%%% - mostly likely culprits are the boundary conditions and boundary cells

%% simulation parameters
% spatial grid
Lx = 100; 
Ly = 100;
nx = 100;
ny = 100;
dx = Lx/nx;
dy = Ly/ny;
x = dx:dx:Lx;
y = dy:dy:Ly;
% time stepping
tmax = 2000;
niter = 20000;
dt = tmax/niter;
t = dt:dt:tmax;

%% physical parameters
k = 1; % thermal conductivity, 2.5 W/mK or 0.3  
rho = 1; % density, 2600 kg/m3 or 2000 
cv = 1; % heat capacity, 2000 J/kg?C or 600
sig = 5.67E-8; % stefan boltzmann constant, 5.67E-8
c = 1; % speed of light 3e8 m/s

Rasteroid = 25; % asteroid radius (Bennu = 493 m, a_distance = 1.10 AU)
% da_dt = ?18.95 � 0.10  (Bennu)
theta = 100; % initial angle
theta = theta/360*2*pi;
w = 0.01; % rad/s, angular velocity, positive = counterclockwise
I = 10; % W/m^2, solar flux on asteroid (1367 W/m2 at Earth)
emission = 1;
Tinit = 0;
%Tinit = (I/4/sig)^(1/4); % initial temp of the asteroid
% set equal to expected average temp to speed up convergence

%% check for numerical instability, terminate if unstable
bx = k/rho/cv * dt/(dx*dx);     % stability parameter
by = k/rho/cv * dt/(dy*dy);     % stability parameter
%bx = 0.334;             % manual input, used for exploring the stability
if ((1-2*bx-2*by)<0)
    disp('WARNING: scheme probably unstable!');
    error(['Multiplicative factor = ',num2str(1-2*bx-2*by)]);    
else
    disp('All fine, scheme should be stable.');
end

%% geometry setup
% all cells are asigned logical variables "is_asteroid" and "is_boundary",
% which are needed for later, especially for boundary conditions
% create distance map from the center, assign all R<R_asteroid to the asteroid
xcenter = round(Lx/2);
ycenter = round(Ly/2);
[xmap,ymap] = meshgrid(x,y);
xmap = xmap-xcenter;
ymap = ymap-ycenter;
Rmap = sqrt(xmap.^2+ymap.^2);
is_asteroid = (Rmap < Rasteroid); % circular asteroid

% boundary = cells belonging to the asteroid but next to a non-asteroid cell
is_boundary = zeros(ny,nx,'logical');
is_boundary_angle = zeros(ny,nx);
boundary_indices = zeros(round(2*pi*Rasteroid/dx),2);
boundary_angle = zeros(round(2*pi*Rasteroid/dx),1);
N_bound = 0;
for i=1:nx
    for j=1:ny
        if is_asteroid(j,i)
            if (~is_asteroid(j+1,i) || ~is_asteroid(j-1,i) ||...
                ~is_asteroid(j,i+1) || ~is_asteroid(j,i-1) )
            N_bound = N_bound+1;
            boundary_indices(N_bound,:) = [j,i];
            boundary_angle(N_bound) = abs(atan2(ymap(j,i)/Rasteroid,xmap(j,i)/Rasteroid)) ...
                        + 2*(pi-abs(atan2(ymap(j,i)/Rasteroid,xmap(j,i)/Rasteroid)))*(asin(ymap(j,i)/Rasteroid)<0);
            is_boundary(j,i) = 1;
            is_boundary_angle(j,i) = boundary_angle(N_bound);
            %acos(xmap(j,i)/Rasteroid) + asin(ymap(j,i)/Rasteroid)...
                    %  + 2*(pi-acos(xmap(j,i)/Rasteroid)).*(asin(ymap(j,i)/Rasteroid)<0)
                    %  + 2*(pi-acos(xmap(j,i)/Rasteroid)).*(asin(ymap(j,i)/Rasteroid)<0)
            end
        end
    end
end
% surf(single(is_boundary)); view(0,90); %shading interp;
surf(is_boundary_angle); view(0,90)
boundary_angle = boundary_angle(1:N_bound);
boundary_indices = boundary_indices(1:N_bound,:);
[boundary_angle, sort_indices] = sort(boundary_angle);
boundary_indices(:,1) = boundary_indices(sort_indices,1);
boundary_indices(:,2) = boundary_indices(sort_indices,2);
%boundary_angle_diff = zeros(N_bound,1);
%for i = 1:N_bound-1
%    boundary_angle_diff(i) = boundary_angle(i+1)-boundary_angle(i);
%end
%boundary_angle_diff(N_bound) = boundary_angle(1)+2*pi-boundary_angle(N_bound);
% OR
boundary_angle_diff = circshift(boundary_angle,1)-circshift(boundary_angle,-1);
boundary_angle_diff(1) = boundary_angle(2)+2*pi-boundary_angle(N_bound);
boundary_angle_diff(N_bound) = boundary_angle(1)+2*pi-boundary_angle(N_bound-1);
boundary_angle_diff = abs(boundary_angle_diff)/2;

% testing the "is_irradiated" variable and output
xmap_rot = cos(theta)*xmap + sin(theta)*ymap;
is_irradiated = is_boundary .* (xmap_rot/Rasteroid);
is_irradiated = is_irradiated > 0; % if on the boundary and on the
% right, cell irradiated by the Sun
surf(is_asteroid+single(is_boundary)+single(is_irradiated)); view(0,90);
hold on;
%plot irradiation vector (with correct angle etc)
quiver( xcenter/dx + cos(theta)*(round(Lx/2)/dx-1), ycenter/dy+sin(theta)*(round(Ly/2)/dy-1),...
        -cos(theta)*round(Lx/8)/dx, -sin(theta)*round(Ly/8)/dy, 0,'y','linewidth',2,'MaxHeadSize',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Loop - diffuse, heat, rotate %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rho*cv*dT/dt = k * laplacian(T);
    % Tnew = T + dt*(k/rho/cv)*(Tleft + Tright - 2*Tcenter)/dx^2
    %          + dt*(k/rho/cv)*(Tup   +  Tdown - 2*Tcenter)/dy^2;
    % Tnew(0,0) = (1-2bx-2by)*T(0,0) + bx*(T(-1,0)+T(+1,0))+ by*(T(0,-1)+T(0,+1));
% initial conditions
To = zeros(ny,nx,'single');
%To(is_asteroid==1 & xmap<0) = Tinit;
To(is_asteroid==1) = Tinit;
T = zeros(ny,nx,niter,'single');
T(:,:,1) = To;
Tmid = To; % dummy variable
Ein = zeros(size(t));
Eout = zeros(size(t));
Eth = zeros(size(t));
angle_irrad = zeros(niter,1);
angle_T = zeros(niter,1);

%%% stepwise physics = rotate + irradiate + emit radiation + diffuse (perhaps different splitting is better?)
for iter = 1:niter % time loop    
    
    Ein_tmp = 0;
    Eout_tmp = 0;
    % calculate which cells are irradiated right now
    xmap_rot = cos(theta)*xmap + sin(theta)*ymap;
    ymap_rot = -sin(theta)*xmap + cos(theta)*ymap;
    is_irradiated = is_boundary .* (xmap_rot/Rasteroid);
    is_irradiated = is_irradiated > 0; % if on the boundary and on the
    % right, cell irradiated by the Sun
    
    % boundary: irradiation and emission
    Tmid = T(:,:,iter); % dummy variable, mid-step temperature
    % irradiation
    for it = 1:N_bound
        i = boundary_indices(it,2);
        j = boundary_indices(it,1);
        if is_irradiated(j,i)
                % Intensity scaled for relative incident angle and 
                % real. vs numerical area factor
                Irel = I * cos(theta-boundary_angle(it))*Rasteroid*boundary_angle_diff(it);
                Ix = abs( Irel * cos(boundary_angle(it)) / dy );
                Iy = abs( Irel * sin(boundary_angle(it)) / dx );
                Ein_tmp = Ein_tmp + Irel;
                
                % dT/dx boundary
                if is_asteroid(j,i+1) == 0 % right boundary
                    %Tmid(j,i) = Tmid(j,i-1) + Ix/k*dx;
                    Tmid(j,i) = Tmid(j,i-1) + Ix/k*dx;
                elseif is_asteroid(j,i-1) == 0 % left boundary
                    %Tmid(j,i) = Tmid(j,i+1) + Ix/k*dx;
                    Tmid(j,i) = Tmid(j,i+1) + Ix/k*dx;
                end
                % dT/dy boundary
                if is_asteroid(j+1,i) == 0 % top boundary
                    %Tmid(j,i) = Tmid(j-1,i) + Iy/k*dy;
                    Tmid(j,i) = Tmid(j-1,i) + Iy/k*dy;
                elseif is_asteroid(j-1,i) == 0 % bottom boundary
                    %Tmid(j,i) = Tmid(j+1,i) + Iy/k*dy;
                    Tmid(j,i) = Tmid(j+1,i) + Iy/k*dy;
                end
        end
    end
    Ein(iter) = Ein_tmp;
                
    % emission
    Tmid2 = Tmid; % dummy variable, mid-step temperature
    if emission
    for it = 1:N_bound
        i = boundary_indices(it,1);
        j = boundary_indices(it,2);
                Iout = sig*T(j,i,iter).^4 * Rasteroid*boundary_angle_diff(it);
                Ix = abs ( Iout * cos(boundary_angle(it))^2 / dy );
                Iy = abs ( Iout * sin(boundary_angle(it))^2 / dx );
                Eout_tmp = Eout_tmp + Iout;
                
                % dT/dx boundary
                if is_asteroid(j,i+1) == 0 % right boundary
                    Tmid(j,i) = Tmid(j,i-1) - Ix/k*dx;
                elseif is_asteroid(j,i-1) == 0 % left boundary
                    Tmid(j,i) = Tmid(j,i+1) - Ix/k*dx;
                end
                % dT/dy boundary
                if is_asteroid(j+1,i) == 0 % top boundary
                    Tmid(j,i) = Tmid(j-1,i) - Iy/k*dy;
                elseif is_asteroid(j-1,i) == 0 % bottom boundary
                    Tmid(j,i) = Tmid(j+1,i) - Iy/k*dy;
                end
    end
    end
    Eout(iter) = Eout_tmp;
    Tmid2 = Tmid;
    Tnew = Tmid2; % dummy variable, "mid-step" temperature
    
    % diffuse
    %%% MERGE INTERNAL AND EXTERNAL NODES?
    for i = 1:nx
        for j = 1:ny
            % INTERNAL NODES
            if is_asteroid(j,i) && ~is_boundary(j,i)
                Tnew(j,i) = (1-2*bx-2*by) .* Tmid2(j,i) +...
                       by .* (Tmid2(j-1,i) + Tmid2(j+1,i)) +...
                       bx .* (Tmid2(j,i+1) + Tmid2(j,i-1));
            end
            % EDGE NODES
            if is_asteroid(j,i) && is_boundary(j,i)
                oneMneigh = 1 - bx*(is_asteroid(j,i-1)+is_asteroid(j,i+1)) ...
                              - by*(is_asteroid(j-1,i)+is_asteroid(j+1,i));
                Tnew(j,i) =  oneMneigh * Tmid2(j,i) +...
                       by * is_asteroid(j-1,i) * Tmid2(j-1,i) +...
                       by * is_asteroid(j+1,i) * Tmid2(j+1,i) +...
                       bx * is_asteroid(j,i+1) * Tmid2(j,i+1) +...
                       bx * is_asteroid(j,i-1) * Tmid2(j,i-1);
            end
        end
    end
    T(:,:,iter+1) = Tnew;
    
    Eth(iter) = (dx*dy)*rho*cv*sum(sum(T(:,:,iter+1) - T(:,:,iter)))/dt;
    
    % hot-cold line angle vs. irradiation angle
    angle_irrad(iter) = theta; % angle of irradiation
    % find max/min temp - can vary randomly, so take average position of
    % cells with temperature in the range [(1-e)*Tmax, Tmax] (same for min)
    Tmax = max(max(Tnew(is_asteroid)));
    Tmin = min(min(Tnew(is_asteroid)));
    dT = Tmax-Tmin;
    e = 0.1;
    [rowmax,colmax] = find(Tnew>=Tmax-e*dT & is_asteroid==1); % indices of Tmax cell
    [rowmin,colmin] = find(Tnew<=Tmin+e*dT & is_asteroid==1); % indices of Tmin cell
    jmax = round(mean(rowmax));
    imax = round(mean(colmax));
    jmin = round(mean(rowmin));
    imin = round(mean(colmin));
    ang = atan2((jmax-jmin),(imax-imin)); % angle of Tmax-Tmin line
    angle_T(iter) = ang .* (ang>=0) + (ang + 2*pi) .* (ang<0);% +...
        %2*pi*floor(theta/(2*pi));
    
    % update angle, next iter
    theta = theta + w*dt;
    
end % end of time loop

%%%%%%%%%%%%%%%%%%%%%%%
%%% post processing %%%
%%%%%%%%%%%%%%%%%%%%%%%
%iter_plot = [1, round(niter/100),round(niter/10), round(niter/1)];
%iter_plot = [1, round(0.8*niter),round(0.85*niter),round(0.9*niter), round(0.95*niter),round(niter/1)];
iter_plot = 6*round(niter/10):round(niter/10):niter;
%figure(iter)
%    surf(T(:,:,iter_plot(1)));
%    axis tight manual
%    set(gca,'nextplot','replacechildren');
%    v = VideoWriter('Bennu.avi');
%    open(v);
    
for iter = 1:length(iter_plot)
    figure(iter)
    surf(T(:,:,iter_plot(iter)));
    view(0,90);
    %frame = getframe(gcf);
    title('Temperature field at t = s')
    xlabel('x (m)')
    ylabel('y (m')
        
    % add irradiation and max-min T vectors
    hold on;
    theta = angle_irrad(iter_plot(iter));
    quiver( xcenter/dx + cos(theta)*(round(Lx/2)/dx-1), ycenter/dy+sin(theta)*(round(Ly/2)/dy-1),...
        -cos(theta)*round(Lx/8)/dx, -sin(theta)*round(Ly/8)/dy, 0,'y','linewidth',2,'MaxHeadSize',2);
    thetaT = angle_T(iter_plot(iter));
    quiver( xcenter/dx + cos(thetaT)*(round(Lx/2)/dx-1), ycenter/dy+sin(thetaT)*(round(Ly/2)/dy-1),...
        -cos(thetaT)*round(Lx/8)/dx, -sin(thetaT)*round(Ly/8)/dy, 0,'b','linewidth',2,'MaxHeadSize',2);
    
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    %hold off
end
%close(v);

%%
figure(7)
    plot(T(round(ny/2),:,1));
    hold on;
    plot(T(round(ny/2),:,round(niter/100)));
    hold on;
    plot(T(round(ny/2),:,round(niter/10)));
    hold on;
    plot(T(round(ny/2),:,niter));
    
figure(8)
   plot(t,Ein,'r');
   hold on;
   plot(t,Eout,'b');
   Eth(1) = 0;
   plot(t,Eth,'g');
   Tavr = mean(Tnew(is_asteroid==1));
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Angle irrad and angle tmax-tmin plot %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
angle_irrad = angle_irrad(1:niter);
angle_T = angle_T(1:niter);
plot(t,angle_irrad);
hold on;
plot(t,angle_T);

figure(8)
plot(t,mod(angle_irrad,2*pi));
hold on;
plot(t,angle_T);

figure(9)
da = mod(angle_irrad,2*pi)-angle_T;
plot(t,da,'b')
da = da .* (da>=0) + (da+2*pi) .* (da<0);
da = da/2/pi*360;
%hold on;
%plot(t,da,'r')
plot(t(round(niter/10):niter),da(round(niter/10):niter));
title('Angle(�) vs. time')
xlabel('time')
ylabel('Angle(�)')

%%    
Tnew = T(:,:,round(niter/1));
Tmax = max(max(Tnew(is_asteroid)));
Tmin = min(min(Tnew(is_asteroid)));
dT = Tmax-Tmin;
e = 0.1;
figure(1)
    surf(is_asteroid + single(is_irradiated) + ...
        2*single(Tnew>=Tmax-e*dT) - 2*single(Tnew<=Tmin+e*dT));
    view(0,90);

figure(2)
    [rowmax, colmax] = find(Tnew>=Tmax-e*dT & is_asteroid==1); % indices of Tmax cell
    [rowmin, colmin] = find(Tnew<=Tmin+e*dT & is_asteroid==1); % indices of Tmin cell
    is_minmax = zeros(ny,nx);
    jmax = round(mean(rowmax));
    imax = round(mean(colmax));
    jmin = round(mean(rowmin));
    imin = round(mean(colmin));
    is_minmax(jmax,imax) = 1;
    is_minmax(jmin,imin) = -1;  
    
    surf(is_asteroid  +5*is_minmax + 3*single(is_irradiated))%
            %+ 4*single(Tnew==Tmax)-4*single(Tnew==Tmin));
    view(0,90);
    

%% alternative
Tavr = mean(Tnew(is_asteroid==1));
imax = 0;
jmax = 0;
imin = 0;
jmin = 0;
dTmax = 0;
dTmin = 0;
for i=1:nx
    for j=1:ny
        if is_asteroid(j,i)
            d = sqrt((x(i)-round(Lx/2))^2 + (y(j)-round(Ly/2))^2);
            if Tnew(j,i) > Tavr
                imax = imax + (Tnew(j,i)-Tavr) * i;
                jmax = jmax + (Tnew(j,i)-Tavr) * j;
                dTmax = dTmax + (Tnew(j,i)-Tavr);
            else
                imin = imin + (Tnew(j,i)-Tavr) * i;
                jmin = jmin + (Tnew(j,i)-Tavr) * j;
                dTmin = dTmin + (Tnew(j,i)-Tavr);
            end
        end
    end
end
imax = round(imax/dTmax);
jmax = round(jmax/dTmax);
imin = round(imin/dTmin);
jmin = round(jmin/dTmin);
is_minmax = zeros(ny,nx);
is_minmax(jmax,imax) = 1;
is_minmax(jmin,imin) = -1;

surf(is_asteroid  +5*is_minmax + 3*single(is_irradiated) +...
    4*single(Tnew==Tmax)-4*single(Tnew==Tmin)); view(0,90);
hold on;

%theta = -atan2((jmax-jmin),(imax-imin)); % angle of Tmax-Tmin line
%theta = theta.* (theta>=0) + (theta + 2*pi) .* (theta<0);% +...
    %2*pi*floor(theta/(2*pi));
%quiver( xcenter+cos(theta)*(round(Lx/2)-1), ycenter-sin(theta)*(round(Ly/2)-1),...
 %       -cos(theta)*round(Lx/8), +sin(theta)*round(Ly/8), 0,'b','linewidth',2,'MaxHeadSize',2);