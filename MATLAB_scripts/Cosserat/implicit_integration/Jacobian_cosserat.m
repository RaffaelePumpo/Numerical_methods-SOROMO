function J = Jacobian_cosserat(time, q, dot_q, ddot_q, Config, Const, R)
%% Jacobian_cosserat - Compute the Jacobian of the residual of the Cosserat rod
% This function calculates the Jacobian matrix of the residual for a
% Cosserat rod. The Jacobian matrix represents the sensitivity of the
% residual to variations in the state variables (positions, velocities, and
% accelerations).
%
% Inputs:
%   - time: Current time
%   - q: Vector of positions
%   - dot_q: Vector of velocities
%   - ddot_q: Vector of accelerations
%   - Config: Configuration parameters for the rod
%   - Const: Constants related to the rod
%   - R: Residual vector for the current state
%
% Output:
%   - J: Jacobian matrix

[Kee, Dee] = computeGeneralisedStiffnessDampingMatrices(Const, Config);
a = Config.a;
b = Config.b;

% Choose the numerical computation method:
%   method = 1: Analytical computation
%   method = 2: Numerical forward computation
%   method = 3: Numerical central computation
method = 1; 
 

%% Analytical Computation
if method == 1
    for i = 1:length(q)
        % Reset Vectors
        Delta_q = zeros(length(q), 1);
        Delta_dot_q = zeros(length(q), 1);
        Delta_ddot_q = zeros(length(q), 1);

        % Propagate the variation
        Delta_q(i) = 1;
        Delta_dot_q(i) = a;
        Delta_ddot_q(i) = b;

        % Compute the numerical variation
        [Delta_Lambda_X0, Delta_Qa_X0] = TIDM(time, q, dot_q, ddot_q, Delta_q, Delta_dot_q, Delta_ddot_q, Config, Const);

        % Calculate the Jacobian columns
        J(:, i) = Kee * Delta_q + Dee * Delta_dot_q - Delta_Qa_X0;
    end
%% Numerical computation
else 
    delta = Config.delta;
    n = length(q);
    for i = 1:n
        delta_q = q;
        delta_dot_q = dot_q;
        delta_ddot_q = ddot_q;

        delta_q(i) = q(i) + delta;
        delta_dot_q(i) = dot_q(i) + a * delta;
        delta_ddot_q(i) = ddot_q(i) + b * delta;
        
        % Calculate the forward difference
        if method == 2
            Rd(:, i) = Residual_cosserat(time, delta_q, delta_dot_q, delta_ddot_q, Const, Config);
            J(:, i) = (Rd(:, i) - R) / delta;
        end
        
        % Calculate the central difference
        if method == 3
            delta_q(i) = q(i) - delta;
            delta_dot_q(i) = dot_q(i) - a * delta;
            delta_ddot_q(i) = ddot_q(i) - b * delta;
            
            Rd2(:, i) = Residual_cosserat(time, delta_q, delta_dot_q, delta_ddot_q, Const, Config);
            
            J(:, i) = (Rd(:, i) - Rd2(:, i)) / (2 * delta);
        end
    end
end
end
