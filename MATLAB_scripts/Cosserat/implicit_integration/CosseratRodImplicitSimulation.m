function [t_stack, states_stack] = CosseratRodImplicitSimulation(q, dot_q, ddot_q, Const,Config)

%Store in Config
t_end = Config.t_end;
dt = Config.dt;
Config.h = 0.01;

%   Threshold for the linear solver
r_min = Config.r_min;


%   Prepare stacks for data
t = 0;
t_stack = [t];
states_stack = [q', dot_q', ddot_q'];


%   Time loop
while t<t_end
    %%  Implicit time integration
    [m, n] = size(states_stack);  % Get the size of the matrix
    qin = states_stack(m,1:length(q))';
    dot_qin = states_stack(m,(length(q)+1):(2*length(q)))';
    ddot_qin = states_stack(m,(2*length(q)+1):n)';

    [q,dot_q,ddot_q] = integrator_cosserat(t,qin,dot_qin,ddot_qin,Config,Const,r_min);    
    
    %   Stack current state
    t_stack = [t_stack; t];
    states_stack = [states_stack;
                    q', dot_q', ddot_q'];
    t = t+dt

end


end