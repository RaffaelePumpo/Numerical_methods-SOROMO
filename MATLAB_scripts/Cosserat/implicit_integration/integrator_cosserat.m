function [qk, dot_qk, ddot_qk] = integrator_cosserat(time, qn, dot_qn, ddot_qn, Config, Const, r_min)
%% integrator_cosserat - Cosserat rod state integrator using the Newmark-beta method with Newton-Raphson correction
%
% This function integrates the state of a Cosserat rod using the Newmark-beta
% method with Newton-Raphson correction. It iteratively corrects the
% predicted state to satisfy the residual equation.
%
% Inputs:
%   - time: Current time
%   - qn: Vector of positions at the current time step
%   - dot_qn: Vector of velocities at the current time step
%   - ddot_qn: Vector of accelerations at the current time step
%   - Config: Configuration parameters for the rod
%   - Const: Constants related to the rod
%   - r_min: Minimum residual norm for convergence
%
% Outputs:
%   - qk: Vector of corrected positions
%   - dot_qk: Vector of corrected velocities
%   - ddot_qk: Vector of corrected accelerations

%% Extracting parameters
dt = Config.dt;
a = Config.a;
b = Config.b;
Beta = Config.Beta;
Gamma = Config.Gamma;

%% Prediction
qnp = qn + dt * dot_qn + (1/2 - Beta) * (dt^2) * ddot_qn;
dot_qnp = dot_qn + (1 - Gamma) * dt * ddot_qn;
ddot_qnp = ddot_qn;
R = Residual_cosserat(time, qnp, dot_qnp, ddot_qnp, Const, Config);

%% Correction
qk = qnp;
dot_qk = dot_qnp;
ddot_qk = ddot_qnp;

while (norm(R) > r_min)
    %% Computation Jacobian
    J = Jacobian_cosserat(time, qk, dot_qk, ddot_qk, Config, Const, R);

    %% Update
    q_up = -pinv(J) * R;
    qk = q_up + qk;
    dot_qk = dot_qk + Gamma / (Beta * dt) * q_up;
    ddot_qk = ddot_qk + 1 / (Beta * dt^2) * q_up;
    R = Residual_cosserat(time, qk, dot_qk, ddot_qk, Const, Config);
end
end
