function [dot_y] = DDM (time, x, Config, Const)
%% DDM - Cosserat rod, Joint velocities and accelerations knowing positions and velocities
% This function computes joint velocities and accelerations of a Cosserat rod
% given positions and velocities using the Direct Dynamic
% Model (DDM).

% Extract positions and velocities from the state vector
q = x(1:Const.dim_base);
dot_q = x((Const.dim_base+1):end);
dot_q0 = zeros(length(dot_q),1);
ddot_q0 = zeros(length(dot_q),1);

% Compute IDM for initial conditions and velocities
[Lambda_X0, Q_c] = IDM(time, q, dot_q0, ddot_q0, Config, Const);
[Lambda_X0, Q_v1] = IDM(time, q, dot_q, ddot_q0, Config, Const);
Q_v = Q_v1 - Q_c;

% Compute the mass matrix M and its columns
for i = 1:1:length(q)
    delta_ddot_q = zeros(length(dot_q),1);
    delta_ddot_q(i) = 1;
    [Lambda_X0, M_1] = IDM(time, q, dot_q, delta_ddot_q, Config, Const);
    M(:,i) = M_1 - Q_c;
end

% Calculate the determinant of M
det_M = det(M);

% Display the determinant of M
disp(['Determinant of M = ' num2str(det_M)]);

% Solve for joint accelerations using inverse of M
ddot_q = -inv(M)*(Q_v + Q_c);

% Construct the output vector
dot_y = [dot_q; ddot_q];
end
