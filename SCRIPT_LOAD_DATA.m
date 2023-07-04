%% Load text files into Workspace
single_cell_table = readtable('01_single_cell.txt');
c0 = table2array(single_cell_table(:,1));
Position = table2array(single_cell_table(:,2:4));
Volume = table2array(single_cell_table(:,5));
M1 = table2array(single_cell_table(:,6));
M2 = table2array(single_cell_table(:,7));
M3 = table2array(single_cell_table(:,8));
M4 = table2array(single_cell_table(:,9));
M5 = table2array(single_cell_table(:,10));
Vel = table2array(single_cell_table(:,11));
M1_grad = table2array(single_cell_table(:,12:14));

neighbor_cell_table = readtable('02_neighbor_cell.txt');
neighbor_cell = table2array(neighbor_cell_table);

inlet_massflow_table = readtable('03_inlet_massflow.txt');
inlet_massflow = table2array(inlet_massflow_table);

outlet_massflow_table = readtable('04_outlet_massflow.txt');
outlet_massflow = table2array(outlet_massflow_table);

% Prepared Data. In this case, M1 is scaled with the hydrodynamic residence
% time value of the CFD

rho=1000;   % density of the fluid in kg/m^3

mass_flow = sum(outlet_massflow(:,3));
volume_flow = mass_flow/rho;
V_ges = sum(Volume);
tau_CFD = V_ges/volume_flow;
c0_M1 = [c0 M1/tau_CFD];
c0_c1_flow = neighbor_cell(:,1:3);
c0_x_y_z_V = [c0 Position Volume];
c0_inlet_flow = [inlet_massflow(:,1:3)];
c0_outlet_flow = [outlet_massflow(:,1:3)];

save('DATA_RAW.mat', 'c0', 'Position', 'Volume', 'M1', 'M2', 'M3', 'M4', 'M5', 'Vel', 'M1_grad', 'tau_CFD', 'volume_flow', 'V_ges', 'neighbor_cell','inlet_massflow','outlet_massflow');

save('DATA_PREPARED.mat','c0_M1','c0_c1_flow','c0_x_y_z_V','c0_inlet_flow','c0_outlet_flow','tau_CFD');