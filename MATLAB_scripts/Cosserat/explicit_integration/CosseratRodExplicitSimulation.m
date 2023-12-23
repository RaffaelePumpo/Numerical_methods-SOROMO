function [t_stack, y_stack] = CosseratRodExplicitSimulation(q, dot_q, Const,Config)

%   Prepare stacks for data
dt = Config.dt;
t_end = Config.t_end;
time_span = 0:dt:t_end;

% Initial conditions
y0 = [q; dot_q];
t_stack = zeros(size(time_span));
y_stack = zeros(length(time_span), length(y0));
%   Explicit time integration
[t_stack,y_stack] = ode45(@(t,x) DDM(t,x,Config,Const),time_span, y0);

end