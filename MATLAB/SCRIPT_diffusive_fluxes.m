%% Approximated diffusive fluxes over the faces
% a diffusive velocity based on u_diff = D / distance is defined. The flux
% A*u_diff is implemented in both direcetions and depends only on geometry,
% not on the velocity field so that it is "valid" for all Re-numbers

clear all;

Re_string = 'RE_1';

single_cell_table = readtable(strcat('01_single_cell_',Re_string,'.txt'));
% for the parfor loop: slice the variables into independent vectors, since
% it increases the calculation speed x_y_z --> x, y, z and c0_c1 --> c0,
% c1. Keep the variables, which are transferred to every worker in the
% parfor loop as small as possible. Volume V is not necessary for the loop,
% the cell number is the shifted index of the array entry, e.g. c0 = 0
% equals index 1 in x_y_z or rather x,y,z
x = table2array(single_cell_table(:,2));
y = table2array(single_cell_table(:,3));
z = table2array(single_cell_table(:,4));

neighbor_cell_table = readtable(strcat('02_neighbor_cell_',Re_string,'.txt'));
neighbor_cell = table2array(neighbor_cell_table);

c0 = neighbor_cell(:,1);    % cell ID 0
c1 = neighbor_cell(:,2);    % cell ID 1
% flow = neighbor_cell(:,3)/hdi_properties(80);
A = neighbor_cell(:,4);     % area of face between c0 and c1

n_connectivity_internal = length(c0);

distance = zeros(n_connectivity_internal,1);
% loop over all internal faces. Get the index of both cells c0 and c1 of
% via the row index i. Since the cell enumeration starts with 0, all cell
% numbers must be increased by 1 for the correct index in the coordinate
% vectors
parfor i = 1:n_connectivity_internal
    distance(i) = sqrt((x(c0(i)+1) - x(c1(i)+1))^2 + (y(c0(i)+1) - y(c1(i)+1))^2 + (z(c0(i)+1) - z(c1(i)+1))^2);
end

% molecular diffusion coefficient in m^2/s
D = 10^-9;
diffusive_flux = A.*D./distance;

% save the "geometric" flux. The diffusive flux needs the diffusion
% coefficient, which may be a variable (molecular or effective diffusion
% coefficient). The specified diffusion coeffcient is a variable in the
% function FUNC_VOLUME_FLOW
geometric_flux = A./distance;
save('DATA_GEOMETRIC_FLUX.mat','geometric_flux');