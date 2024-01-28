function [T04,T03] = transformation_T04(th1,th2,th3,th4)
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
% finding the transformation of each coordinate frame wrt to base frame
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T03 = real(double(subs(T03,{theta1,theta2,theta3},{th1,th2,th3})));
T04 = real(double(subs(T04,{theta1,theta2,theta3,theta4},{th1,th2,th3,th4})));