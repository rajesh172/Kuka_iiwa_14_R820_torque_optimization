clc;
clear all;
close all;
load tau_max_without_null_space.mat
Force_ellipsoid=[]; % data for force ellipsoid
dt=1/10;
t=0:dt:2*pi;%timesteps
% defining the trajectory which is infinity in 3d-space
y = 0.2*cos(t);
z = cos(t).*sin(2*t)*0.2+0.078663+0.78;
x = t./(pi*12)-0.526;
% finding the endeffector velocity
dx=diff(x)./dt;
dy=diff(y)./dt;
dz=diff(z)./dt; 
v=[dx;dy;dz];

q=[0;0;0;pi/2;0;0;0];  % initial condition
x_ = [0;0;0]; % for plotting the endeffector coordinates given by forward kinematics
F = [1;1;1];
tau_max =[];
tau =[];

NO = [0.1;0.3;0.2;0.15;-0.2;0.1;-0.5]; % random null space vector

for k=1:1:length(x)-1
    th1=q(1,k);th2=q(2,k);th3=q(3,k);th4=q(4,k);th5=q(5,k);th6=q(6,k);th7=q(7,k);
    Jv = Jacobian_force_6_7(th1,th2,th3,th4,th5,th6,th7);  
    tau(:,k) = Jv'*F;
    tau_max(k) = sqrt(tau(1,k)^2+ tau(2,k)^2+ tau(3,k)^2+tau(4,k)^2+tau(5,k)^2+tau(6,k)^2+tau(7,k)^2);
    [xo,yo,zo] = fwd_kin_kuka_iiwa_14R820(q(:,k));
    x_(:,k)=[xo;yo;zo];

    E=[x(k);y(k);z(k)]-[xo;yo;zo];% error btw actual trajectory coordinates and coordinates from forward kinematics    
    % Updating the angle 
    
    % NC = (eye(7)-pinv(Jv)*Jv)*NO; % null-space projection


    q(:,k+1)=q(:,k) + 1*pinv(Jv)*(v(:,k) + 0.2*E)*dt + dt*(eye(7)-pinv(Jv)*Jv)*NO(:,1);

    % drawing the the force ellipsoid at the three points on trajectory
    % lower value, middle value and higher value of sqrt(det(force_ellips))
    force_ellips = pinv(Jv*Jv');
    temp = sqrt(det(force_ellips));
    Force_ellipsoid = [Force_ellipsoid;temp];
    % k=1 value of sqrt(det(force_ellips)) is low, at k=38 -it is in mid range, and at k=58 -it is high 
    if k ==1 || k==38 || k ==58
       ellips = force_ellips;
       [evc,evl] = eig(ellips);
       figure();
       axis equal;
       ellipsoid(0,0,0,evl(1,1),evl(2,2),evl(3,3));% drawing ellipsoid taking centre at (0,0,0) at particular point on trajectory
       xlim([-15,15]);
       ylim([-15,15]);
       title(["force ellipsoid at k =",num2str(k)])
       xlabel("x-axis");
       ylabel("y-axis")
       zlabel("z-axis");
    end
end
t(end)=[];
figure;
plot(t,tau_max,"Color",[1,0,0])
hold on
plot(t,tau_max_without_null_space',"Color",[0,0,1])
xlabel('time-steps')
ylabel('Torque')
title('Null-space effect')
legend('tau with null-space','tau without null-space')

