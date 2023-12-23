function Const = defineCosseratRod(Config)



%%  GEOMETRICAL PARAMETERS

%   Lenght of the rod
Const.L = 1;

%   Radius Cross section
Const.Rc  = 0.001;

%   Area
Const.Area = pi*Const.Rc^2;

%   Geometrical moment of inertia
Const.J      = zeros(3,3);
Const.J(1,1) = pi*Const.Rc^4/2;
Const.J(2,2) = pi*Const.Rc^4/4;
Const.J(3,3) = pi*Const.Rc^4/4;



%%  MATERIAL PARAMETERS

%   Specific weight
Const.rho = 7800;

%   Shear modulus
Const.G = 80e9;

%   Young modulus
Const.E  = 210e9;

%   Internal damping
Const.mu = 1e-4;



%%  STIFFNESS AND INERTIA OF CROSS SECTION

Const.GIxx = Const.G*Const.J(1,1);
Const.EIyy = Const.E*Const.J(2,2);
Const.EIzz = Const.E*Const.J(3,3);
Const.EA = Const.E*Const.Area;
Const.GA = Const.G*Const.Area;

Const.GI = Const.GIxx;
Const.EI = Const.EIyy;


Const.M      = zeros(3,3);
Const.M(1,1) = Const.rho*pi*Const.Rc^2;
Const.M(2,2) = Const.rho*pi*Const.Rc^2;
Const.M(3,3) = Const.rho*pi*Const.Rc^2;

Const.M_cal  = [Const.rho*[Const.J,zeros(3)];[zeros(3),Const.M]];



%%  STRAIN BASED PARAMETERIZATION


% Definition of actuated values
% K1, K2, K3, ...
% if value is 1 -> the variable is actuated
% if value is 0 -> the variable is not actuated
Const.V_a = [0, 1, 0, 0, 0, 0];

%   Define the size of the parameterisation
ne = 3;
Const.dim_base_k = ne*Const.V_a;
%   Compute size of q
Const.dim_base   = Const.V_a*Const.dim_base_k';


%   Automatically define matrix B
M_selec = eye(6,6);
[~,col] = find(Const.V_a==1);
Const.B = M_selec(:,col);

%   Define constant strain
Const.Xi_c = [0;0;0;1;0;0];



%%  STIFFNESS AND DAMPING MATRICES

[Kee, Dee] = computeGeneralisedStiffnessDampingMatrices(Const, Config);

Const.Kee = Kee;
Const.Dee = Dee;



end