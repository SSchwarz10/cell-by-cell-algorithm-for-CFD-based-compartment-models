function [] = FUNC_MOMENTS_LINEAR_SYSTEM(n_SLICES,n_LOCAL,Diff,Re)
% n_SLICES: number of slices for the domain
% n_LOCAL: divisions of local Mean-Age within a slice
% Diff: Diffusion coefficient for the diffusive flux
% Re: Reynoldsnumber for the operating point

% Reynolds-number for the base case / CFD
Re_base = 1;

% preparing the variable strings for loading, saving, the diary and the
% interpolation file
load_string = strcat('CLUSTERED_DATA_RE_1_VOLUME_FLOW_SLICES_',num2str(n_SLICES),'_DELTA_',num2str(n_LOCAL),'.mat');
save_string = strcat('RESULT_RE_BASE_1_D_',num2str(Diff),'_RE_',num2str(Re),'_MOMENTS.mat');
diary_string = strcat('Diary_RE_BASE_1_D_',num2str(Diff),'_RE_',num2str(Re),'_moments.txt');
interpolation_string = strcat('INTERPOLATION_RE_BASE_1_D_',num2str(Diff),'_RE_',num2str(Re),'_MOMENTS.txt');

tic

% load the data
load(load_string,'c0_x_y_z_V','CPT','CPT_slice','CPT_V',...
    'Flow_Outlet','tau','n_CPT','Flow_Matrix_conv','Flow_Matrix_geometric');

% create a fileID for writing the diary
fileID = fopen(diary_string,'w');

Flow_Matrix_diff = Flow_Matrix_geometric*Diff;

% scaling factor
scaling = Re/Re_base;

% the overall Flow_Matrix is the sum of the convective and the diffusive
% flow matrix

Flow_Matrix = Flow_Matrix_conv*scaling + Flow_Matrix_diff;

% the matrix does not contain any entries on the diagoanl yet, which
% represent the total ouflow of every compartment. Due to the mass balance
% for every compartment and the assumption that the density is constant,
% the overall outgoing flux for a compartment is equal to the sum of all
% ingoing fluxes. This is represented by the sum along every row. Since the
% flow is an outgoing flux, the sign is negative

D = diag(-sum(Flow_Matrix,2));
Flow_Matrix = Flow_Matrix + D;

% alle fluxes are volumetric flow rates. Since the balance for a CSTR is
% based on the volume, every row must be "scaled" with the corresponding
% compartment volume (./Volume_all). The inlet and outlet are associated
% with a volume of 1. Thus, only the last two rows in the matrix contains
% the exact volume flow rates. Additionally, since there are no diffusive
% fluxes from or to the inlet or outelt, the last two rows only contains
% the convective streams

Volume_all = [CPT_V(:,2);1;1];
Flow_Normalized = Flow_Matrix./Volume_all;

% since the transport equations for the moments can be treated as a
% reaction of 0-th order for every moment, the ode is actually a system of
% linear equations: dM_n/dt = Q * M_n + n*M_(n-1) =!= 0
% Thus, the vector of the moments can be obtained by inverting the flow
% matrix Q and conduct a simple matrix calculation. Since the flow matrix
% is the same for every moment, it is useful to invert it once and use it
% up to the moment of interest. Special care must be taken for the
% Flow_Normalized matrix: since the system is closed-closed, there is no
% flow to the "inlet" compartment and no flow from the "outlet" compartment
% into the domain, the before last and last column/row may not be used.
% Otherwise, the determinant of the matrix is 0 and can not be inverted.

Q_rev = inv(Flow_Normalized(1:end-2,1:end-2));

M1_all = -Q_rev * ones(n_CPT,1);
M2_all = -Q_rev * 2 * M1_all;
M3_all = -Q_rev * 3 * M2_all;
M4_all = -Q_rev * 4 * M3_all;
M5_all = -Q_rev * 5 * M4_all;

Flow_Avg_M1 = sum(M1_all.*(Flow_Outlet.'))/sum(Flow_Outlet);
Flow_Avg_M2 = sum(M2_all.*(Flow_Outlet.'))/sum(Flow_Outlet);
Flow_Avg_M3 = sum(M3_all.*(Flow_Outlet.'))/sum(Flow_Outlet);
Flow_Avg_M4 = sum(M4_all.*(Flow_Outlet.'))/sum(Flow_Outlet);
Flow_Avg_M5 = sum(M5_all.*(Flow_Outlet.'))/sum(Flow_Outlet);

% steady state of the moments
M1 = Flow_Avg_M1(end);
M2 = Flow_Avg_M2(end);
M3 = Flow_Avg_M3(end);

M1 = full(M1);
M2 = full(M2);
M3 = full(M3);

% calcualte the caracterisitc key numbers
variance = M2 - M1^2;
sigma2_theta = variance/M1^2;
CoV = sqrt(variance/M1^2);
central_M3 = M3 - 3*M2*M1 + 2*M1^3;
skewness = (central_M3/M1^3)^(1/3);
math_skewness = central_M3/((M2-M1^2)^(3/2));

% calculate the charactersitics for every slice for a comparison with the
% raw values of the CFD
% extract the steady state value for every moment (:) for all compartments
% (1:end-2)
M1_end = M1_all;
M2_end = M2_all;
M3_end = M3_all;
M4_end = M4_all;
M5_end = M5_all;

Slice_M1_M2_M3_M4_M5_var_CoV = zeros(n_SLICES,8);
for k = 1:n_SLICES
    CPT_within_slice = CPT_slice(CPT_slice(:,2)==k,1);
    V_within_slice = CPT_V(CPT_V(CPT_within_slice,1),2);
    M1_within_slice = M1_end(CPT_within_slice,:);
    M2_within_slice = M2_end(CPT_within_slice,:);
    M3_within_slice = M3_end(CPT_within_slice,:);
    M4_within_slice = M4_end(CPT_within_slice,:);
    M5_within_slice = M5_end(CPT_within_slice,:);
    M1_Vavg = sum(M1_within_slice.*V_within_slice)/sum(V_within_slice);
    M2_Vavg = sum(M2_within_slice.*V_within_slice)/sum(V_within_slice);
    M3_Vavg = sum(M3_within_slice.*V_within_slice)/sum(V_within_slice);
    M4_Vavg = sum(M4_within_slice.*V_within_slice)/sum(V_within_slice);
    M5_Vavg = sum(M5_within_slice.*V_within_slice)/sum(V_within_slice);
    
    var_within_slice = sum((M1_within_slice - M1_Vavg).^2.*V_within_slice)/sum(V_within_slice);
    CoV_within_slice = sqrt(var_within_slice)/M1_Vavg;
    
    Slice_M1_M2_M3_M4_M5_var_CoV(k,1) = k;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,2) = M1_Vavg;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,3) = M2_Vavg;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,4) = M3_Vavg;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,5) = M4_Vavg;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,6) = M5_Vavg;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,7) = var_within_slice;
    Slice_M1_M2_M3_M4_M5_var_CoV(k,8) = CoV_within_slice;
end
% Volume average of M1
M1_Vavg = sum(M1_end.*CPT_V(:,2))/sum(CPT_V(:,2));
% Volume average of M2
M2_Vavg = sum(M2_end.*CPT_V(:,2))/sum(CPT_V(:,2));
% Volume average of M1^2
M1_sqr_Vavg = sum(M1_end.^2.*CPT_V(:,2))/sum(CPT_V(:,2));

% CoV regarding age, mean age and degree of mixing Zwietering and
% Danckwerts (defined as ratio of vaiances w.r.t to volume averaged age)
CoV_age_ZD = sqrt(M2_Vavg - M1_Vavg^2)/M1_Vavg;
CoV_mean_age_ZD = sqrt(M1_sqr_Vavg - M1_Vavg^2)/M1_Vavg;
J_ZD = (M1_sqr_Vavg - M1_Vavg^2)/(M2_Vavg - M1_Vavg^2);

% CoV regarding age, mean age and degree of mixing Liu (defined as ratio of
% vaiances w.r.t to mean residence time)

CoV_age_L = sqrt(M2_Vavg - 2*M1_Vavg*tau + tau^2)/tau;
CoV_mean_age_L = sqrt(M1_sqr_Vavg - 2*M1_Vavg*tau + tau^2)/tau;
% J_L = (M1_sqr_Vavg - 2*M1_Vavg*tau + tau^2) / (M2_Vavg - 2*M1_Vavg*tau + tau^2);
J_L = (M1_sqr_Vavg - 2*M1_Vavg*M1 + M1^2) / (M2_Vavg - 2*M1_Vavg*M1 + M1^2);

%%
% writing interpolation file for FLUENT to show the moment values of the
% compartment model in the CFD mesh. Export the moment value as the same
% UDS as in FLUENT. Inputmatrix must have the form : x,y,z,M1,M2,M3,M4,M5
% last time step (end) for all compartments (:). transposed to a column
% vector (.')
M_end = [M1_end M2_end M3_end M4_end M5_end];
% 3 coordinates + n_COMPONENTS
Inputmatrix = zeros(size(CPT,1), 3 + 5);
Inputmatrix(:,1:3) = c0_x_y_z_V(:,2:4);
for j = 1:n_CPT
    % divide the entries with j since CPT(CPT==j) is the compartment number
    % j
    Inputmatrix(CPT==j,4:end) = CPT(CPT==j).*M_end(j,:)./j;
end

species = {'uds-0' 'uds-1' 'uds-2' 'uds-3' 'uds-4'};

exporttofluent(single(Inputmatrix),species,interpolation_string);

time=toc;
toc;

% fprintf(fileID,'Simulated time [s] [s]: %d\n',t(end));
fprintf(fileID,'Time for calculation [s]: %d\n',time);
fclose(fileID);
save(save_string,'Flow_Avg_M1','Flow_Avg_M2','Flow_Avg_M3','Flow_Avg_M4','Flow_Avg_M5',...
    'tau','M_end','M1','M2','M3','variance','sigma2_theta','CoV','central_M3','skewness','math_skewness',...
    'n_SLICES','n_LOCAL','n_CPT','Slice_M1_M2_M3_M4_M5_var_CoV','CPT_slice',...
    'CoV_age_ZD','CoV_mean_age_ZD','J_ZD','CoV_age_L','CoV_mean_age_L','J_L')