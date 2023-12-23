clc
clear all
close all
% Start the timer
tic;
%   Time integration settings
t_end = 10;
dt = 1e-2;
time_span = 0:dt:t_end;
Beta  = 0.25;
Gamma = 0.5;
a = Gamma/(Beta*dt);
b = 1/(Beta*dt^2);

%Store in Config
Config.dt = dt;
Config.a = a;
Config.b = b;
Config.Beta = Beta;
Config.Gamma = Gamma;


%   Threshold for the linear solver
r_min = 1e-6;


%   Properties of the system
Config.m = 1;
Config.l = 1;
Config.k = 10;
Config.g = 9.81;
%   Initial conditions
qr_0 = -pi/4;
qp_0 = 0.5;
% qr_0 = 3;
% qp_0 = 1;
q = [qr_0; qp_0];
dot_q  = [0,0]';
ddot_q = [0,0]';

%   Prepare stacks for data
t = 0;
t_stack = [t];
states_stack = [q', dot_q', ddot_q'];


%   Time loop
while t<t_end

    %%  Implicit time integration
    [m, n] = size(states_stack);  
    y = states_stack(m,:);
    [q,dot_q,ddot_q] = implicit_residual(t,y,Config,r_min);
    
    %   Stack current state
    t_stack = [t_stack; t];
    states_stack = [states_stack;
                    q', dot_q', ddot_q'];
    t = t+dt

end

% Stop the timer
elapsed_time = toc;

% Display the elapsed time
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
%   Plot results
figure("Name", "Joint Values2")
plot(t_stack, states_stack(:,1), '-r', 'DisplayName', 'q1')
hold on 
plot(t_stack, states_stack(:,2), '-b', 'DisplayName', 'q2')
grid on 
axis equal
legend('Location', 'best')



%%
for it_t =1:max(size(t_stack))
    q1 = states_stack(it_t,1);
    q2 = states_stack(it_t,2);

    l = Config.l;

    m_x = (l + q2)*cos(q1);
    m_y = (l + q2)*sin(q1);

    figure(1)
    plot([0 m_x], [0, m_y],'r','linewidth',2)
    hold on

    grid on
    hold off
    axis equal
    xlim([-3 3])
    ylim([-3 3])
    title(['time = ' num2str(t_stack(it_t)) 's'])
    drawnow
end


