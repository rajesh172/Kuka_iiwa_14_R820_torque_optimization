clc;clear all;close all;
sim=remApi('remoteApi'); % using the prototype file (remoteApiProto.m)
sim.simxFinish(-1); % just in case, close all opened connections
clientID = sim.simxStart('127.0.0.1',19999,true,true,5000,5);

if (clientID>-1)
    disp('Connected to remote API server');
    %creates some joint pos
    dt=1/10;
    penLinearPosition_x = []; % endeffector position from copeliasim
    penLinearPosition_y= [];
    penLinearPosition_z = [];
    % joints handles of our manipulator kuka from the copeliasim
    h = [0,0];
    [err,h(1)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/joint',sim.simx_opmode_blocking);
    [err,h(2)]= sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/link2_resp/joint',sim.simx_opmode_blocking);
    [err,tar] = sim.simxGetObjectHandle(clientID,'/LBRiiwa14R820/Pen',sim.simx_opmode_blocking);
    i=1;
    for q1=-170*pi/180:dt:170*pi/180

        for q2=-120*pi/180:dt:120*pi/180
            % Set joint positions of our manipulator in copeliasim
            sim.simxSetJointTargetPosition(clientID,h(1),q1,sim.simx_opmode_streaming);       
            sim.simxSetJointTargetPosition(clientID,h(2),q2,sim.simx_opmode_streaming);
    
            % Get pen(endeffector) linear position
            [err, LinearPosition_ee] = sim.simxGetObjectPosition(clientID, tar,-1, sim.simx_opmode_blocking);
            penLinearPosition_x(i) = LinearPosition_ee(1); % x-coordinates of the position of endeffector
            penLinearPosition_y(i) = LinearPosition_ee(2); % y-coordinates of the position of endeffector
            penLinearPosition_z(i) = LinearPosition_ee(3); % z-coordinates of the position of endeffector
            i=i+1;
            pause(0.05); 
        end
        i=i+1;

    end
    figure;
    plot3(penLinearPosition_x,penLinearPosition_y,penLinearPosition_z);
    title('workspace of kuka iiwa 14R820');
    xlabel("x"); ylabel("y"); zlabel("z");

else
    disp('Failed connecting to remote API server'); % shown when the matlab is not connected with copeliasim
end
sim.delete(); % call the destructor!
disp('Program ended');


