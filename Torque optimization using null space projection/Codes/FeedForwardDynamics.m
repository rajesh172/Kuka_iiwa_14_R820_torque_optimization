% Equation of motion
function [tau1,tau2,tau3,tau4,tau5,tau6,tau7] = FeedForwardDynamics(q,qd,qdd)
% torque for jooint 1
tau1 = subs(T1_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 2
tau2 = subs(T2_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 3
tau3 = subs(T3_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 4
tau4 = subs(T4_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 5
tau5 = subs(T5_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 6
tau6 = subs(T6_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});

% torque for joint 7
tau7 = subs(T7_subs,{t1,t2,t3,t4,t5,t6,t7,t1d,t2d,t3d,t4d,t5d,t6d,t7d,t1dd,t2dd,t3dd,t4dd,t5dd,t6dd,t7dd}, ...
    {q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7)});
