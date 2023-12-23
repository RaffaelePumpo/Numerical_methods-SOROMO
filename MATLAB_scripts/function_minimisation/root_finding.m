function root = root_finding(x0, f, method)
%% root_finding - Find the root of the function f using various methods
%
% This function aims to find the root of the given function f using
% different methods specified by the 'method' parameter.
%
% Inputs:
%   - x0: Initial guess for the root
%   - f: Function handle representing the equation to solve
%   - method: Method for root finding (1 for analytical, 2 for numerical forward, 3 for numerical central)
%
% Output:
%   - root: Approximate root of the function f

% Initializations
xv = x0;
dx = 1e-6;
R = f(xv);
epsilon = 1e-6;
i = 0;

% Iterative root-finding process
while (norm(R) > epsilon && i < 1000)

    %% Analytical
    if method == 1
        J = @(x) 4*x^3 - 12*x^2 + 10*x - 2;
        J = J(xv);
    end

    %% Numerical forward
    if method == 2
        Rd = f(xv+dx);
        J = (Rd - R) / dx;
    end

    %% Numerical central
    if method == 3
        Rd = f(xv+dx);
        Rd2 = f(xv - dx);
        J = (Rd - Rd2) / (2 * dx);
    end

    % Update root estimate
    xup = -pinv(J) * R;
    xv = xv + xup;
    R = f(xv);
    i = i + 1;
end

% Output the approximate root
root = xv;
end
