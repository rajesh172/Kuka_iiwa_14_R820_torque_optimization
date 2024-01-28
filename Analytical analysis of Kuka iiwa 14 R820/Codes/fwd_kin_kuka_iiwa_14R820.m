function [x,y,z] = fwd_kin_kuka_iiwa_14R820(q)
 % lengths of the different links as given in kuka dh parameter
  d1 = 0.36;
  d3 = 0.42;
  d5 = 0.4 ;
  d7 = 0.126;
syms theta alph d a;
% general matrix used when findinf the jacobian from DH parameter
M = [cos(theta), -sin(theta)*cos(alph), sin(theta)*sin(alph), a*cos(theta);sin(theta),cos(theta)*cos(alph),-cos(theta)*sin(alph), a*sin(theta); 0,sin(alph),cos(alph),d; 0, 0, 0, 1];

% substituting the values of DH parameter form each joint in the above
% matrix
A1 = subs(M,{a,alph,d,theta},{0,-pi/2,d1,q(1,:)});
A2 = subs(M,{a,alph,d,theta},{0,pi/2,0,q(2,:)});
A3 = subs(M,{a,alph,d,theta},{0,pi/2,d3,q(3,:)});
A4 = subs(M,{a,alph,d,theta},{0,-pi/2,0,q(4,:)});
A5 = subs(M,{a,alph,d,theta},{0,-pi/2,d5,q(5,:)});
A6 = subs(M,{a,alph,d,theta},{0,pi/2,0,q(6,:)});
A7 = subs(M,{a,alph,d,theta},{0,0,d7,q(7,:)});

% finding the transformation of each coordinate frame wrt to base frame
T01 = A1;
T02 = T01*A2;
T03 = T02*A3;
T04 = T03*A4;
T05 = T04*A5;
T06 = T05*A6;
T07 = T06*A7;% trasformation of endeffector frame wrt base frame

% end-effector coordinates from the first 3 rows of the 4th column of the
% T07 transformation matrix
x = double(T07(1,4));
y = double(T07(2,4));
z = double(T07(3,4));
end