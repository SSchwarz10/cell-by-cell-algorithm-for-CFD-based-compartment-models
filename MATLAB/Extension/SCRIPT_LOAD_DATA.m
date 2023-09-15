%% Load text files into Workspace

Re_string = 'RE_1';

Re = 1;         % reynolds-number in 1
rho=1000;       % density in kg/m^3
dyn_vis = 10^-3;    % dynamic viscosity in Pa*s
kin_vis = dyn_vis/rho;  % kinematic viscosity in m^2/s
d_pipe = 12.4*10^-3;    % pipe diameter in m
u_superficial = Re*kin_vis/d_pipe;    % superficial velocity in m/s
D=10^-9;        % diffusion coefficient in m^2/s
Sc= kin_vis/D;  % schmidt-number in -

single_cell_table = readtable(strcat('01_single_cell_',Re_string,'.txt'));
c0 = table2array(single_cell_table(:,1));
Position = table2array(single_cell_table(:,2:4));
Volume = table2array(single_cell_table(:,5));
M1 = table2array(single_cell_table(:,6));
M2 = table2array(single_cell_table(:,7));
M3 = table2array(single_cell_table(:,8));
M4 = table2array(single_cell_table(:,9));
M5 = table2array(single_cell_table(:,10));
% Variance = table2array(single_cell_table(:,11));
% Rec = table2array(single_cell_table(:,12));
Vel = table2array(single_cell_table(:,11));
M1_grad = table2array(single_cell_table(:,12:14));

neighbor_cell_table = readtable(strcat('02_neighbor_cell_',Re_string,'.txt'));
neighbor_cell = table2array(neighbor_cell_table);

inlet_massflow_table = readtable(strcat('03_inlet_massflow_',Re_string,'.txt'));
inlet_massflow = table2array(inlet_massflow_table);

outlet_massflow_table = readtable(strcat('04_outlet_massflow_',Re_string,'.txt'));
outlet_massflow = table2array(outlet_massflow_table);

% Prepared Data. In this case, M1 is scaled with the hydrodynamic residence
% time value of the CFD
mass_flow = sum(outlet_massflow(:,3));
volume_flow = mass_flow/rho;
V_ges = sum(Volume);
tau_CFD = V_ges/volume_flow;
c0_M1 = [c0 M1/tau_CFD];
c0_c1_flow = neighbor_cell(:,1:3);
c0_x_y_z_V = [c0 Position Volume];
c0_inlet_flow = [inlet_massflow(:,1:3)];
c0_outlet_flow = [outlet_massflow(:,1:3)];

save(strcat('DATA_RAW_',Re_string,'.mat'), 'c0', 'Position', 'Volume',...
    'M1', 'M2', 'M3', 'M4', 'M5', 'Vel', 'M1_grad',...
    'tau_CFD', 'volume_flow', 'V_ges', 'neighbor_cell','inlet_massflow','outlet_massflow',...
    'Re','Sc','u_superficial','d_pipe','rho','dyn_vis','kin_vis','D');

save(strcat('DATA_PREPARED_',Re_string,'.mat'),'c0_M1','c0_c1_flow','c0_x_y_z_V','c0_inlet_flow','c0_outlet_flow','tau_CFD',...
    'Re','Sc','u_superficial','d_pipe','rho','dyn_vis','kin_vis','D');