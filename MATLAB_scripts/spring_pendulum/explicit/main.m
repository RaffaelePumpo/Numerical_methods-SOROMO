clc
clear all
close all
tic;
%   Time integration settings
t_end = 10;
dt = 1e-2;
time_span = 0:dt:t_end;

%   Properties of the system
Config.m = 1;
Config.l = 1;
Config.g = 9.81;
Config.k = 10;

%   Initial conditions
%qr_0 = -pi/4;
%qp_0 = 0.5;
%  qr_0 = pi/2;
% qp_0 = 1;
dot_qr_0 = 0;
dot_qp_0 = 0;
 % dot_qr_0 = 3;
 % dot_qp_0 = 3;
y0 = [qr_0; qp_0; dot_qr_0; dot_qp_0];


%   Explicit time integration
[t,y]=ode45(@(t,x) DDM(t,x,Config),time_span, y0);
elapsed_time = toc;

% Display the elapsed time
disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);
% %   Plot results
figure("Name", "Joint Values")
plot(t, y(:,1), '-r', 'DisplayName', 'q1')
hold on 
plot(t, y(:,2), '-b', 'DisplayName', 'q2')
grid on 
axis equal
legend('Location', 'best')

drawnow


%%
for it_t =1:max(size(t))
    q1 = y(it_t,1);
    q2 = y(it_t,2);

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
    title(['time = ' num2str(t(it_t)) 's'])
    drawnow
end




