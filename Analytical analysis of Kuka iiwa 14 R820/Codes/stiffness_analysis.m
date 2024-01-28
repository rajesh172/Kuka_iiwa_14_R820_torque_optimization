clc;clear all;close all;
stif_ellips=[]; % data for stiffness 

dt=1/10;
t=0:dt:2*pi; %timesteps
numSteps = length(t);

% defining the trajectory which is infinity in 3d-space
y = 0.2*cos(t);
z = cos(t).*sin(2*t)*0.2+0.078663+0.78;
x = t./(pi*12)-0.526;
% defining the joint space stiffness matrix
kq =[0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25;
    0.1,0.2,0.05,0.6,0.4,0.15,0.25];
% defining the linear force on the endffector
F = [1;0;0];

% finding the endeffector velocity
dx=diff(x)./dt;
dy=diff(y)./dt;
dz=diff(z)./dt; 

v=[dx;dy;dz]; % velocity vector
q=[0;0;0;pi/2;0;0;0];   % initial configuration
x_ = [0;0;0]; % for plotting the endeffector coordinates given by forward kinematics

for k=1:1:length(x)-1
    th1=q(1,k);th2=q(2,k);th3=q(3,k);th4=q(4,k);th5=q(5,k);th6=q(6,k);th7=q(7,k);
    [Jv,kg] = Jacobian_kuka_with_stiff_kg(th1,th2,th3,th4,th5,th6,th7,F);% jacobian
    [xo,yo,zo] = fwd_kin_kuka_iiwa_14R820(q(:,k)); % finding endeffector coordinates from forward kinematics
    x_(:,k)=[xo;yo;zo];
    E=[x(k);y(k);z(k)]-[xo;yo;zo]; % error btw actual trajectory coordinates and coordinates from forward kinematics    
    % Updating the angle 
    q(:,k+1)=q(:,k) + 1*pinv(Jv)*(v(:,k) + 0.2*E)*dt;

    % calculating the task space stiffness matrix using CCT(conservative
    % congruency transformation)
    kx = pinv(Jv')*(kq-kg)*pinv(Jv);

    % drawing the the stiffness ellipsoid at the three points on trajectory
    % lower value, middle value and higher value of sqrt(det(stiffness))
    stiffness = pinv(kx*kx');
    temp = sqrt(det(stiffness));
    stif_ellips = [stif_ellips;temp];
    % k=56 value of sqrt(det(force_ellips)) is low, at k=5 -it is in mid range, and at k=30 -it is high 
    if k ==5 || k==56 || k ==30
       manipula = stiffness;
       [evc,evl] = eig(manipula);
       figure();
       axis equal;
       ellipsoid(0,0,0,evl(1,1),evl(2,2),evl(3,3)); % drawing ellipsoid taking centre at (0,0,0) at particular point on trajectory
       xlim([-4,4]);
       ylim([-4,4]);
       title(["stiffness ellipsoid at k =",num2str(k)])
       xlabel("x-axis");
       ylabel("y-axis")
       zlabel("z-axis");
    end
end
% since our stiff_ellips vector is lesser length because of
% differentiation, that's why removing last element from x,y,z
x(end)=[];
y(end)=[];
z(end)=[];
%% 

figure();
% using triscatteredInterp function and surf to plot the heat map in three
% dimension to have better understandin at which point on trajectory the value of
% sqrt(det(stiffness)) is higher or lower
F = TriScatteredInterp(y',z',stif_ellips);
ti = -0.3:.01:0.9;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
s = surf(qx,qy,qz);
s.EdgeColor = 'none';
axis square ;
% axis equal ;
grid on ;
hold on
title('','Interpreter','Latex') ;
xlabel('x-coordinate','FontSize',20,'Interpreter','Latex') ;
ylabel('y-coordinate','FontSize',20,'Interpreter','Latex') ;
zlabel('Stiffness-index','FontSize',20,'Interpreter','Latex') ;
set(gca,'fontsize',20) ;
set(gca,'LineWidth',1) ;

