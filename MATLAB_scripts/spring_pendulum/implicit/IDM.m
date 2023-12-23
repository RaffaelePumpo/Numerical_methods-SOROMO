function tau = IDM(t,y,Config)
%% More stable wrt Explicit that is sensible to stiff systems that tends to diverce 
%% This is more elaborate but more stable
%% impose tau = taud st 
%% Based on prediction and correctiono
%% Residual is tau-taud = 0
q1 = y(1);
q2 = y(2);
dot_q1 = y(3);
dot_q2 = y(4);
ddot_q =[y(5),y(6)]';


m = Config.m;
l = Config.l;
g = Config.g;
k = Config.k;
M = [m*(l+q2)^2, 0 ;
        0,         m];

C = [2*m*dot_q1*dot_q2*(l+q2);
    -m*dot_q1^2*(l+q2)];
K = [m*g*(l+q2)*cos(q1);
    m*g*sin(q1)+k*q2];

tau = M*ddot_q + C+K;
end