function [qk,dot_qk,ddot_qk] = implicit_residual(t,y,Config,r_min)
method = 1;
%% Extracting parameters
dt = Config.dt;
Beta = Config.Beta;
Gamma = Config.Gamma;

qn = y(1:2)';

dot_qn = y(3:4)';

ddot_qn = y(5:6)';

%% Prediction
qnp = qn+dt*dot_qn + (1/2 -Beta)*(dt^2)*ddot_qn;
dot_qnp = dot_qn + (1-Gamma)*dt*ddot_qn;
ddot_qnp = 0*ddot_qn;
y = [qnp',dot_qnp',ddot_qnp'];
R = IDM(t,y,Config);

%% Correction
qk = qnp;
dot_qk = dot_qnp;
ddot_qk = ddot_qnp;
while (norm(R)>r_min)

    %% Computation Jacobian
    J = Jacobian_spring(qk,dot_qk,ddot_qk,Config,R,method);

    %% Update
    q_up = -pinv(J)*R;
    qk = q_up + qk;
    dot_qk = dot_qk + Gamma/(Beta*dt)*q_up;
    ddot_qk = ddot_qk + 1/(Beta*dt^2)*q_up;
    y = [qk',dot_qk',ddot_qk'];
    R = IDM(t,y,Config);
end
end