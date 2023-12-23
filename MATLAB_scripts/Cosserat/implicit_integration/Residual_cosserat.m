function R = Residual_cosserat(time, q, dot_q, ddot_q, Const, Config)
%% Residual_cosserat - Compute the residual for a Cosserat rod
% This function calculates the residual for a Cosserat rod based on the
% provided positions, velocities, accelerations, and system parameters.

% Compute the generalised stiffness and damping matrices
[Kee, Dee] = computeGeneralisedStiffnessDampingMatrices(Const, Config);

% Compute the internal damping and mass matrices using IDM
[Lambda_X0, Qa_X0] = IDM(time, q, dot_q, ddot_q, Config, Const);

% Calculate the residual using the stiffness, damping, and internal forces
R = Kee * q + Dee * dot_q - Qa_X0;
end
