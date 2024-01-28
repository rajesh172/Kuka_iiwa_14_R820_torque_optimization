clc;clear all;close all;
sim=remApi('remoteApi'); % using the prototype file (remoteApiProto.m)
sim.simxFinish(-1); % just in case, close all opened connections
clientID = sim.simxStart('127.0.0.1',19999,true,true,5000,5);
manipulability_ellips=[];
if (clientID>-1)
    disp('Connected to remote API server');
    %creates some joint pos
    dt=1/10;
    t=0:dt:2*pi; %timesteps
    numSteps = length(t);
    penLinearVelocities = zeros(3, numSteps); %endeffector velocity  from copeliasim
    penLinearPosition = zeros(3, numSteps); % endeffector position from copeliasim
    % defining the trajectory which is infinity in 3d-space
    y = 0.2*cos(t);
    z = cos(t).*sin(2*t)*0.2+0.078663+0.78;
    x = t./(pi*12)-0.526;

    % finding the endeffector velocity
    dx=diff(x)./dt;
    dy=diff(y)./dt;
    dz=diff(z)./dt; 
    v=[dx;dy;dz]; % velocity vector
    q=[0;0;0;pi/2;0;0;0];   % initial configuration
    x_ = [0;0;0]; % for plotting the endeffector coordinates given by forward kinematics

    for k=1:1:length(x)-1
        th1=q(1,k);th2=q(2,k);th3=q(3,k);th4=q(4,k);th5=q(5,k);th6=q(6,k);th7=q(7,k);
        Jv = Jacobian_kuka_iiwa_14R820(th1,th2,th3,th4,th5,th6,th7);% jacobian
        [xo,yo,zo] = fwd_kin_kuka_iiwa_14R820(q(:,k)); % finding endeffector coordinates from forward kinematics
        x_(:,k)=[xo;yo;zo];
        E=[x(k);y(k);z(k)]-[xo;yo;zo]; % error btw actual trajectory coordinates and coordinates from forward kinematics  

        
        % Updating the angle 
        q(:,k+1)=q(:,k) + 1*pinv(Jv)*(v(:,k) + 0.2*E)*dt;
        
        % drawing the the manipulability ellipsoid at the three points on trajectory
        % lower value, middle value and higher value of sqrt(det(force_ellips))
        manipulability = Jv*Jv';
        temp = sqrt(det(manipulability));
        manipulability_ellips = [manipulability_ellips;temp];
        % k=1 manipulability is high, k=38, it is in mid range and at k=58 manipulability is low  
        if k ==1 || k==38 || k ==58
           manipula = manipulability;
           [evc,evl] = eig(manipula);
           figure();
           axis equal;
           ellipsoid(0,0,0,evl(1,1),evl(2,2),evl(3,3));% drawing ellipsoid taking centre at (0,0,0) at particular point on trajectory
           xlim([-4,4]);
           ylim([-4,4]);
           title(["Manipulability ellipsoid at k =",num2str(k)])
           xlabel("x-axis");
           ylabel("y-axis")
           zlabel("z-axis");
        end
    end
    figure;
    plot3(x_(1,:),x_(2,:),x_(3,:)); %plotting the endeffector coordinates given by forward kinematics
    title('tajectory from forward and inverse kinematics');
    xlabel("x"); ylabel("y"); zlabel("z");
  

    % joints handles of our manipulator kuka from the copeliasim
    h = [0,0,0,0,0,0,0];
    [err,h(1)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/joint',sim.simx_opmode_blocking);
    [err,h(2)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link2_resp/joint',sim.simx_opmode_blocking);
    [err,h(3)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link3_resp/joint',sim.simx_opmode_blocking);
    [err,h(4)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link4_resp/joint',sim.simx_opmode_blocking);
    [err,h(5)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link5_resp/joint',sim.simx_opmode_blocking);
    [err,h(6)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link6_resp/joint',sim.simx_opmode_blocking);
    [err,h(7)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link7_resp/joint',sim.simx_opmode_blocking);
    [err,tar] = sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/Pen',sim.simx_opmode_blocking);

    for i=1:numSteps
        % Set joint positions of our manipulator in copeliasim
        sim.simxSetJointTargetPosition(clientID,h(1),q(1,i),sim.simx_opmode_streaming);       
        sim.simxSetJointTargetPosition(clientID,h(2),q(2,i),sim.simx_opmode_streaming);
        sim.simxSetJointTargetPosition(clientID,h(3),q(3,i),sim.simx_opmode_streaming);
        sim.simxSetJointTargetPosition(clientID,h(4),q(4,i),sim.simx_opmode_streaming);
        sim.simxSetJointTargetPosition(clientID,h(5),q(5,i),sim.simx_opmode_streaming);
        sim.simxSetJointTargetPosition(clientID,h(6),q(6,i),sim.simx_opmode_streaming);
        sim.simxSetJointTargetPosition(clientID,h(7),q(7,i),sim.simx_opmode_streaming);

        % Get pen(endeffector) linear velocity
        [err, LinearVelocity_ee] = sim.simxGetObjectVelocity(clientID, tar, sim.simx_opmode_blocking);
        penLinearVelocities(1,i) = LinearVelocity_ee(1); % x-coordinates of the velocity of endeffector
        penLinearVelocities(2,i) = LinearVelocity_ee(2); % y-coordinates of the velocity of endeffector
        penLinearVelocities(3,i) = LinearVelocity_ee(3); % z-coordinates of the velocity of endeffector

        % Get pen(endeffector) linear position
        [err, LinearPosition_ee] = sim.simxGetObjectPosition(clientID, tar,-1, sim.simx_opmode_blocking);
        penLinearPosition(1,i) = LinearPosition_ee(1); % x-coordinates of the position of endeffector
        penLinearPosition(2,i) = LinearPosition_ee(2); % y-coordinates of the position of endeffector
        penLinearPosition(3,i) = LinearPosition_ee(3); % z-coordinates of the position of endeffector

        %pause(0.05);  

    end

    % Get simulation time
    simTime = t; % Simulation time is synchronized with the loop

    figure;
    % Plot pen linear position of endeffector
    subplot(2, 1, 1);
    plot(simTime, penLinearPosition);
    title('Pen Linear Position vs. Time');
    xlabel('Time (s)');
    ylabel('Position');
    legend('Position_x','Positions_y','Position_z');
    % % Plot pen linear velocity
    subplot(2, 1, 2);
    plot(simTime, penLinearVelocities);
    title('Pen Linear Velocity vs. Time');
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    legend('velocity_x','velocity_y','velocity_z');

    % Plot joint angle of the manipulator wrt time just to see whether we
    % are getting correct angles from copeliasim
    figure;
    plot(simTime, q);
    title('joint angle vs. Time');
    xlabel('Time (s)');
    ylabel('joint angle');
    legend('joint1','joint2','joint3','joint4','joint5','joint6','joint7');
    % plotting the endeffector actual position that is trace by kuka in the
    % copeliasim
    figure;
    plot3(penLinearPosition(1,:),penLinearPosition(2,:),penLinearPosition(3,:));
    title('tajectory from copeliasim');
    xlabel("x"); ylabel("y"); zlabel("z");

else
    disp('Failed connecting to remote API server'); % shown when the matlab is not connected with copeliasim
end
sim.delete(); % call the destructor!
disp('Program ended');

%%
% using triscatteredInterp function and surf to plot the heat map in three
% dimension to have better understandin at which point on trajectory the value of
% sqrt(det(stiffness)) is higher or lower
x(end)=[];
y(end)=[];
z(end)=[];
figure();
F = TriScatteredInterp(y',z',manipulability_ellips);
ti = -0.3:.01:0.9;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
s = surf(qx,qy,qz);
s.EdgeColor = 'none';
axis square ;
% axis equal ;
grid on ;
hold on
title('manipulability index','Interpreter','Latex') ;
xlabel('y-coordinate','FontSize',20,'Interpreter','Latex') ;
ylabel('z-coordinate','FontSize',20,'Interpreter','Latex') ;
zlabel('Manipulability-index','FontSize',20,'Interpreter','Latex') ;
set(gca,'fontsize',20) ;
set(gca,'LineWidth',1) ;

