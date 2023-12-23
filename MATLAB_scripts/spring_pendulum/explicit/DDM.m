function dot_y = DDM(t,y,Config)
%% DDM Pendolum, joint accelerations and velocities knowing  position and velocity
q1 = y(1);
q2 = y(2);
dot_q1 = y(3);
dot_q2 = y(4);
m = Config.m;
l = Config.l;
g = Config.g;
k = Config.k;
%% Considering no force and no torque
tau = [0,0]';

%% Model
M = [m*(l+q2)^2, 0 ;
        0,         m];

C = [2*m*dot_q1*dot_q2*(l+q2);
    -m*dot_q1^2*(l+q2)];
K = [m*g*(l+q2)*cos(q1);
    m*g*sin(q1)+k*q2];
ddot_q = inv(M)*(tau-C-K);

dot_y = [dot_q1;
         dot_q2;
         ddot_q];
end