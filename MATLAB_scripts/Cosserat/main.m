close all
clear all
clc


addpath("ODEs/")
addpath("utilities/")
addpath("rod_properties/")
addpath("implicit_integration/")
addpath("explicit_integration/")


%   Configurations for the simulation
Config = simulationConfigurations();

%   Define properties for the rod
Const = defineCosseratRod(Config);

%   Value of gravity
Const.g = 9.81;

%   Position of the rod base
Const.r_X0 = [0;0;0];

%   Quaternion of rod base [w x y z]
%Const.Q_X0 = [0.7071068 0 0.7071068 0]';
Const.Q_X0 = [1, 0, 0, 0]';

%   Cosserat rod generalized coordinates
q      = zeros(Const.dim_base, 1);
dot_q  = zeros(Const.dim_base, 1);
ddot_q = zeros(Const.dim_base, 1);


%% Implicit Computation
% Start the timer
tic;
[t_stack_implicit, states_stack_implicit] = CosseratRodImplicitSimulation(q, dot_q, ddot_q, Const,Config);
% Stop the timer
elapsed_time = toc;
% Display the elapsed time
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);

%% Explicit Computation
tic;
[t_stack_explicit, states_stack_explicit] =  CosseratRodExplicitSimulation(q, dot_q, Const,Config); 
% Stop the timer
elapsed_time = toc;
% Display the elapsed time
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);

%% Plot joint 

% figure("Name", "Joint Values")
% plot(t_stack_implicit, states_stack_implicit(:,1), '-r', 'DisplayName', 'q1')
% hold on 
% plot(t_stack_implicit, states_stack_implicit(:,2), '-b', 'DisplayName', 'q2')
% hold on 
% plot(t_stack_implicit, states_stack_implicit(:,3), '-k', 'DisplayName', 'q3')
% grid on 
% axis equal
% legend('Location', 'best')

%% Plot implicit

 for it_t=1:1   %length(t_stack_implicit)

    t = t_stack_implicit(it_t);

    ne = Const.dim_base;

    q = states_stack_implicit(it_t, 1:ne)';

    plotRod(t, q, Const, Config)

end

%% Plot Explicit

for it_t=1:max(size(t_stack_explicit))

    t = t_stack_explicit(it_t);

    ne = Const.dim_base;

    q = states_stack_explicit(it_t, 1:ne)';

    plotRod(t, q, Const, Config)

end



