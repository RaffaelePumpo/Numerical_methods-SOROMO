function Config = simulationConfigurations()

%   Time integration settings
t_end = 3;
dt    = 1e-2;
Beta  = 1/4;
Gamma = 1/2;
a     = Gamma/(Beta*dt);
b     = 1/(Beta*dt^2);


%   Store in Config
Config.t_end = t_end;
Config.dt    = dt;
Config.a     = a;
Config.b     = b;
Config.Beta  = Beta;
Config.Gamma = Gamma;

%   Properties for Cosserat ODEs integration
Config.dX = 10^-3;
Config.forward_integration_domain = 0:Config.dX:1;
Config.backward_integration_domain = 1:-Config.dX:0;


%   Threshold residual norm
Config.r_min = 1*10^-7;


%   Value of the numerical perturbation
Config.delta = 1e-6;

end