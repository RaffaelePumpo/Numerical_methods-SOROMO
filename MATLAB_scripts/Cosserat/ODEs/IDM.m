function [Lambda_X0, Qa_X0] = IDM(time, q, dot_q, ddot_q, Config, Const)


%   Define the forward state
forward_y0 = zeros(19, 1);

%   Initialization
forward_y0(1:7) = [Const.Q_X0;
                    Const.r_X0];

%   Forward integration
[~, Y] = ode45(@(X, y) ForwardKinematics(X, y, q, dot_q, ddot_q, Config, Const), ...
                    Config.forward_integration_domain, forward_y0);

%   Extract state at X=1
y_tip = Y(end,:);

Q_X1       = y_tip(1:4)';
r_X1       = y_tip(5:7)';
eta_X1     = y_tip(8:13)';
dot_eta_X1 = y_tip(14:19)';

%   Force at the tip
F1 = F_tip(time);

%   Define the forward state
backward_y1 = zeros(26 + Const.dim_base - 1, 1);

%   Initialization
backward_y1(1:19) = [Q_X1;
                     r_X1;
                     eta_X1;
                     dot_eta_X1];
backward_y1(20:25) = F1;

%   Backward integration
[~, Y] = ode45(@(X, y) BackwardDynamics(X, y, q, dot_q, ddot_q, Config, Const), ...
                    Config.backward_integration_domain, backward_y1);

%   Extract state at X=0
y_base = Y(end,:);

Lambda_X0 = y_base(20:25)';
Qa_X0 = y_base(26:end)';


end

