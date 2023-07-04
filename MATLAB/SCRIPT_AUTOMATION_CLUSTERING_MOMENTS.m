%% Script for creating and calculating of the different compartment models for one case

n_SLICES = [38];   % number of slices for the domain
n_LOCAL = [1];    % divisions of local Mean-Age within a slice

rho =   1000;   % density of the fluid in kg/m^3

% "DATA_PREPARED"  must be in the same directory like this script to create
% the compartment model, determine the flow matrix and apply the mean-age
% theory to the CPT model

tic
% calling the functions within 
for i = 1:length(n_SLICES)
    for j = 1:length(n_LOCAL)
        % clustering the CFD-cells
        FUNC_CLUSTERING(n_SLICES(i),n_LOCAL(j));
        % calculating the exchange streams
        FUNC_VOLUME_FLOW(n_SLICES(i),n_LOCAL(j),rho);
        % calculating the mean age and higer moments field
        FUNC_MOMENTS_LINEAR_SYSTEM(n_SLICES(i),n_LOCAL(j));
    end
end
toc

% Outcome:
% CLUSTERED_DATA_SLICES_XX_DELTA_YY:
% contains the data for the clustered CFD cells
% CLUSTERED_DATA_VOLUME_FLOW_SLICES_XX_DELTA_YY:
% contains the flow matrix for the CPT model and the TIS model
% RESULT_SLICES_38_DELTA_1_MOMENTS:
% contains the mean-age and its higher moments field and the
% characteristics of the global RTD based on the mean age