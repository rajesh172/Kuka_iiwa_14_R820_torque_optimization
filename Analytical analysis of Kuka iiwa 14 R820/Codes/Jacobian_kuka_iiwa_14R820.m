function J = Jacobian_kuka_iiwa_14R820(th1,th2,th3,th4,th5,th6,th7)
  % lengths of the different links as given in kuka dh parameter
  d1 = 0.36;
  d3 = 0.42;
  d5 = 0.4 ;
  d7 = 0.126;

syms theta alph d a theta1 theta2 theta3 theta4 theta5 theta6 theta7;
% general matrix used when findinf the jacobian from DH parameter
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];

% substituting the values of DH parameter form each joint in the above
% matrix
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,d1,theta1});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta2});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,d3,theta3});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,theta4});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,d5,theta5});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,theta6});
A7 = subs(M,{a,alph,d,theta},{0,0,d7,theta7});

% finding the transformation of each coordinate frame wrt to base frame
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06*A7; % trasformation of endeffector frame wrt base frame

% defining the jacobian by differentiating the end-effector coordinates wrt
% to each joing angle
Jv = [diff(T07(1,4),theta1),diff(T07(1,4),theta2),diff(T07(1,4),theta3),diff(T07(1,4),theta4),diff(T07(1,4),theta5),diff(T07(1,4),theta6),diff(T07(1,4),theta7);...
    diff(T07(2,4),theta1),diff(T07(2,4),theta2),diff(T07(2,4),theta3),diff(T07(2,4),theta4),diff(T07(2,4),theta5),diff(T07(2,4),theta6),diff(T07(2,4),theta7);...
    diff(T07(3,4),theta1),diff(T07(3,4),theta2),diff(T07(3,4),theta3),diff(T07(3,4),theta4),diff(T07(3,4),theta5),diff(T07(3,4),theta6),diff(T07(3,4),theta7)];

% replacing the numerical values of thetas in jacobian
J = double(subs(Jv,{theta1,theta2,theta3,theta4,theta5,theta6,theta7},{th1,th2,th3,th4,th5,th6,th7}));
end
