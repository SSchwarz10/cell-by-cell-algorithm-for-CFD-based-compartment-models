%% Script for creating and calculating of the different compartment models for one case

n_SLICES = [42];   % number of slices for the domain
n_LOCAL = [15];    % divisions of local Mean-Age within a slice

Diff = 1*10^-9;       % diffusion coefficient in m^2/s

rho = 1000;    % density in kg/m^3
tic
% calling the functions within 
for i = 1:length(n_SLICES)
    for j = 1:length(n_LOCAL)
        % clustering the CFD-cells
        FUNC_CLUSTERING(n_SLICES(i),n_LOCAL(j));
        % calculating the exchange streams
        FUNC_VOLUME_FLOW_CONV_DIFF(n_SLICES(i),n_LOCAL(j),rho);
    end
end

n_Re = [1 2 5 10 20 50 100 200];

for i = 1:length(n_Re)
    % calculate the moments for the different operating points
    FUNC_MOMENTS_LINEAR_SYSTEM(n_SLICES(1),n_LOCAL(1),Diff,n_Re(i));
end


toc