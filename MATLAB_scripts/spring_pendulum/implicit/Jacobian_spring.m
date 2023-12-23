function J = Jacobian_spring(qk,dot_qk,ddot_q_k,Config,R,method)
%Function to compute Jacobian of the spring pendolum
% method =1 is the analytical computation of the jacobian
% method =2 is the  forward numerical computation of the jacobian 
% method =1 is the central numerical computation of the jacobian
    a = Config.a;
    b = Config.b;
%% Analytically
if method == 1
    m = Config.m;
    g = Config.g;
    l= Config.l;
    k= Config.k;
    qk1 = qk(1);
    qk2 = qk(2);
    dot_qk1 = dot_qk(1);
    dot_qk2 = dot_qk(2);
    ddot_qk1 = ddot_q_k(1);

    dtau1q1 = m*b*((l+qk2)^2)+2*a*m*dot_qk2*(l+qk2)-m*g*(l+qk2)*sin(qk1);
    dtau1q2 = 2*m*ddot_qk1*(l+qk2)+2*a*m*dot_qk1*(l+qk2)+2*m*dot_qk1*dot_qk2+m*g*cos(qk1);
    dtau2q1 = -2*a*m*(l+qk2)*dot_qk1+m*g*cos(qk1);
    dtau2q2 = m*b-m*(dot_qk1^2) + k;
    J = [dtau1q1, dtau1q2;
    dtau2q1,dtau2q2];
%% Numerically
else 
    delta = 1e-6;
    for i = 1:1:2
        delta_q = qk;
        delta_dot_q = dot_qk;
        delta_ddot_q = ddot_q_k;

        delta_q(i) =  qk(i) + delta;
        delta_dot_q(i) = dot_qk(i) + a*delta;
        delta_ddot_q(i) = ddot_q_k(i) + b*delta;
        y = [delta_q',delta_dot_q',delta_ddot_q'];
        Rd(:,i) = IDM(1,y,Config);
        %% Forward
        if method == 2
            J(:,i) = (Rd(:,i)-R)/delta;
        end
        %% Central
        if method == 3
        delta_q = qk;
        delta_dot_q =dot_qk;
        delta_ddot_q = ddot_q_k;

        delta_q(i) =  qk(i) - delta;
        delta_dot_q(i) = dot_qk(i) - a*delta;
        delta_ddot_q(i) = ddot_q_k(i) - b*delta;
        y = [delta_q',delta_dot_q',delta_ddot_q'];
        Rd2(:,i) = IDM(1,y,Config);
        J(:,i) = (Rd(:,i) - Rd2(:,i))/(2*delta);
        end


    end

end
end