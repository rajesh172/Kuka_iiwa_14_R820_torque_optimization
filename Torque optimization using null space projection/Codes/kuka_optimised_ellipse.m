clc;
%clear;
clc
%%
clearvars -except opti
for ci=6:1:6
load('NR_x.mat')

q=[0,0,0,-pi/2,0,0,0];
dt=1/100;
t=0:dt:2*pi;

Datai = double(forward_kinematics_kuka(q));
%x = t./(pi*12);
x=ones(size(t))*Datai(1,1);
y = cos(t)*0.2;
z = sin(t)*0.2;
%x=x+Datai(1,1) - x(1,1);
y=y+Datai(1,2)-y(1,1);
z=z+Datai(1,3)-z(1,1);
dx= diff(x)/dt;
dy = diff(y)/dt;
dz = diff(z)/dt;
v=[dx;dy;dz;dx*0;dy*0;dz*0];
Av=[1E-6,1E-6,1E-6,1E-6,1E-6,1E-6];
Av(ci)=-1;
q=q';
clear Data
Data=zeros(1,6);
% b_i=0.04;
% c_i=0.04;
% d_i=0.04;
eff=[];
for k=1:1:length(x)-1
    th1=q(1,k);th2=q(2,k);th3=q(3,k);th4=q(4,k);th5=q(5,k);th6=q(6,k);th7=q(7,k);
    J =double(Jacobian_kuka_iiwa(q(:,k))); 
    disp(k)
    Do = double(forward_kinematics_kuka(q(:,k)'));
    Data(k,:)=Do;
    % show(lbr, q(:,k)', 'PreservePlot', false);
    % hold on
    % plot3(Data(1:k,1),Data(1:k,2),Data(1:k,3),'r','linewidth',1)
    % [X,Y,Z]=TDellipse(J(1:3,:),Do);
    % draw_circle(Do,ci);
    % quiver3(Do(1),Do(2),Do(3),Av(1)/2,Av(2)/2,Av(3)/2,'linewidth',1.5,'color','b');
    
    axis([-1.5 1.5 -1.5 1.5 0 2])
    hold off
    drawnow
    set(gca)
    F(k) = getframe(gcf);
    
    E=[x(k);y(k);z(k);0;pi/2;0]-Do';
    
    no = nullspce_cotribution_qo(th1,th2,th3,th4,th5,th6,th7,Av(1),Av(2),Av(3),Av(4),Av(5),Av(6));
    NC=(eye(7) - pinv(J)*J)*no ;
    q(:,k+1)=q(:,k) + 1*pinv(J)*(v(:,k) + 1E2*E)*dt + 0*1E2*NC*dt;
    tau=J'*[Av(1);Av(2);Av(3);Av(4);Av(5);Av(6)];
    eff(k,:)=sqrt(tau'*tau);


    
    
end
%% 

%   writerObj = VideoWriter(strcat('optimised_circle',num2str(ci),'.avi'));
%   %writerObj = VideoWriter(strcat('without_optimised_circle','.avi'));
%   writerObj.FrameRate = 30;
% %   writerObj.Width = 560;
% % writerObj.Height = 420;
% 
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ; 
% 
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);
end 
%%

function []= draw_circle(Do,ci)
t=0:0.01:pi;
if ci==4
    x = Do(1)*ones(size(t));
    y = ones(size(t))*Do(2) + 0.1*sin(t);
    z = ones(size(t))*Do(3) + 0.1*cos(t);
    plot3([x(end) x(end)],[y(end) y(end)+0.04],[z(end) z(end)-0.04],'b','linewidth',1.5)
    plot3([x(end) x(end)],[y(end) y(end)+0.04],[z(end) z(end)+0.04],'b','linewidth',1.5)
    plot3(x,y,z,'b','linewidth',1.5);
end
if ci==5
    y = Do(2)*ones(size(t));
    x = ones(size(t))*Do(1) + 0.1*sin(t);
    z = ones(size(t))*Do(3) + 0.1*cos(t);
    plot3([x(end) x(end)+0.04],[y(end) y(end)],[z(end) z(end)+0.04],'b','linewidth',1.5)
    plot3([x(end) x(end)+0.04],[y(end) y(end)],[z(end) z(end)-0.04],'b','linewidth',1.5)
    plot3(x,y,z,'b','linewidth',1.5);
end
if ci==6
    z = Do(3)*ones(size(t));
    x = ones(size(t))*Do(1) + 0.1*sin(t);
    y = ones(size(t))*Do(2) + 0.1*cos(t);
    plot3([x(end) x(end)+0.04],[y(end) y(end)+0.04],[z(end) z(end)],'b','linewidth',1.5)
    plot3([x(end) x(end)+0.04],[y(end) y(end)-0.04],[z(end) z(end)],'b','linewidth',1.5)
    plot3(x,y,z,'b','linewidth',1.5);
end



end

function [X,Y,Z]=TDellipse(J,Do)

% Calculate the force ellipsoid matrix
velocity_ellipsoid_matrix = J*J';

% Perform eigenvalue decomposition
[V, D] = eig(velocity_ellipsoid_matrix);
D=(inv(D)./max(max(inv(D))))./4;

% Generate points on the unit sphere
[u, v] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi/2, pi/2, 50));
x = cos(u) .* cos(v);
y = sin(u) .* cos(v);
z = sin(v);

% Scale and rotate the points on the unit sphere
unit_sphere_points = [x(:), y(:), z(:)]';
velocity_ellipsoid_points = V * sqrt(D)* unit_sphere_points;

% Reshape the transformed points for plotting
X = reshape(velocity_ellipsoid_points(1, :), size(x));
Y = reshape(velocity_ellipsoid_points(2, :), size(y));
Z = reshape(velocity_ellipsoid_points(3, :), size(z));
X = X + ones(50,50)*Do(1);
Y = Y + ones(50,50)*Do(2);
Z = Z + ones(50,50)*Do(3);
surf(X, Y, Z,'EdgeColor','none','FaceAlpha',0.25);
end

function [Data] = forward_kinematics_kuka(q)
syms theta alph d a;
q=q';
a1 = 0.36;
a2 = 0.42;
a3 = 0.4;
a4 = 0.126;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,a1,q(1,:)});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,q(2,:)});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,a2,q(3,:)});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,q(4,:)});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,a3,q(5,:)});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,q(6,:)});
A7 = subs(M,{a,alph,d,theta},{0,0,a4,q(7,:)});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06*A7;
R = [T07(1,1),T07(1,2),T07(1,3);T07(2,1),T07(2,2),T07(2,3);T07(3,1),T07(3,2),T07(3,3)];
x = T07(1,4);
y = T07(2,4);
z = T07(3,4);
angle_y = pi/2;
angle_x = 0;
angle_z = 0;
Data = [x,y,z,angle_x,angle_y,angle_z];
end


function J = Jacobian_kuka_iiwa(q)
syms theta alph d a theta1 theta2 theta3 theta4 theta5 theta6 theta7;
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,0.36,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,0.42,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,0.4,theta5});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta6});
A7 = subs(M,{a,alph,d,theta},{0,0,0.126,theta7});
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06*A7;

z0 = [0; 0; 1];
z1 = [T01(1,3);T01(2,3);T01(3,3)];
z2 = [T02(1,3);T02(2,3);T02(3,3)];
z3 = [T03(1,3);T03(2,3);T03(3,3)];
z4 = [T04(1,3);T04(2,3);T04(3,3)];
z5 = [T05(1,3);T05(2,3);T05(3,3)];
z6 = [T06(1,3);T06(2,3);T06(3,3)];


Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4),diff(T07(1,4),theta5),diff(T07(1,4),theta6),diff(T07(1,4),theta7);
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4),diff(T07(2,4),theta5),diff(T07(2,4),theta6),diff(T07(2,4),theta7);
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4),diff(T07(3,4),theta5),diff(T07(3,4),theta6),diff(T07(3,4),theta7)];
Jw = [z0 z1 z2 z3 z4 z5 z6];
J = [Jv;Jw];
J = subs(J,{theta1,theta2,theta3,theta4,theta5,theta6,theta7},{q(1),q(2),q(3),q(4),q(5),q(6),q(7)});
end

