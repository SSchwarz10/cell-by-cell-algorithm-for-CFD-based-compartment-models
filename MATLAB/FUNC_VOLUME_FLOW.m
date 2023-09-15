function [] = FUNC_VOLUME_FLOW(n_SLICES,n_LOCAL,rho)
%% function for the overall clustering process, which is called within the script "SCRIPT_AUTOMATION_CLUSTERING_MOMENTS". The data for load are created by the function "FUNC_CLUSTERING"

% n_SLICES: number of slices for the domain
% n_LOCAL: divisions of local Mean-Age within a slice
% rho: density for calculating the volume flow from the mass flow

% preparing the variable strings for loading, saving and the diary
load_string = strcat('CLUSTERED_DATA_SLICES_',num2str(n_SLICES),'_DELTA_',num2str(n_LOCAL),'.mat');
save_string = strcat('CLUSTERED_DATA_VOLUME_FLOW_SLICES_',num2str(n_SLICES),'_DELTA_',num2str(n_LOCAL),'.mat');
diary_string = strcat('Diary_slices_',num2str(n_SLICES),'_delta_',num2str(n_LOCAL),'_flowrates.txt');

% load the data
load(load_string,'c0_M1','c0_c1_flow','c0_inlet_flow','c0_outlet_flow','CPT_V','c0_x_y_z_V','Barycenter','Vol_Avg_M1','CPT','delta_x','n_SLICES','n_LOCAL','CPT_slice','n_CPT','n_CELLS');

% create a fileID for writing the diary
fileID = fopen(diary_string,'w');
fprintf(fileID,('\n-----START: Calculation of exchange streams-----\n'));
tic
%%
% copy c0_c1_flow to overwrite the cell number with the compartment number.
% Similar procedure for the inlet and outlet matrices
CPT_flow = c0_c1_flow;
n_connectivity_internal = size(CPT_flow,1);

for i = 1:n_connectivity_internal
    for j =1:2
        if CPT_flow(i,j) ~= -1 && CPT_flow(i,j) ~=-2
            CPT_flow(i,j) = CPT(CPT_flow(i,j)+1);
        end
    end
end

CPT1_inlet_flow = c0_inlet_flow;
n_connectivity_inlet = size(CPT1_inlet_flow,1);

for i = 1:n_connectivity_inlet
    CPT1_inlet_flow(i,1) = CPT(CPT1_inlet_flow(i,1)+1);
end

CPT1_outlet = c0_outlet_flow;
n_connectivity_outlet = size(CPT1_outlet,1);

for i = 1:n_connectivity_outlet
    CPT1_outlet(i,1) = CPT(CPT1_outlet(i,1)+1);
end
%%
% matrix for storing the flow from CPT1 to CPT2 and vice versa. Entries
% with a zero in the first or second column will be deleted
CPT1_CPT2_flow=zeros(n_connectivity_internal,3);

%fluxes within one compartment are irrelevant and are deleted to reduce
%loop iterations
CPT_flow(CPT_flow(:,1) == CPT_flow(:,2),:) = [];

% counter for the rows in CPT1_CPT2_flow
row = 1;

% copy of CPT_flow for manipulation
CPT_flow_copy = CPT_flow;

% end of time measurement
t=toc;

fprintf(fileID,'Time for preparation [s]: %d\n',t);

% loop over all compartnemts
for i = 1:n_CPT
    start = toc;
    % all rows containing the current compartment number are copied
    CPT_flow_needed = CPT_flow_copy(:,1)==i | CPT_flow_copy(:,2)==i;
    CPT_flow_small = CPT_flow_copy(CPT_flow_needed,:);
    
    % identification of all neighbored compartments
    CPT_unique = unique(CPT_flow_small(:,1:2));
    
    % delelte the current compartment number in order to reduce the loop
    % iterations and get the number of neighbored compartments
    CPT_unique(CPT_unique == i) = [];
    n_neighbor = size(CPT_unique,1);
    % loop over all neighbors
    for j = 1:n_neighbor
        % identification of all required rows
        % neighbored compartment B
        B = CPT_unique(j);
        % there is a flux from compartment i to B if these compartments are
        % in the same order within a row an the flux in the third column is
        % positive OR the order of i and B are switched and the flux in the
        % third row is negative. This condition is met in one logical array
        % array_i_B
        array_i_B = (CPT_flow_small(:,1)==i & CPT_flow_small(:,2)==B & CPT_flow_small(:,3)>=0) | (CPT_flow_small(:,2)==i & CPT_flow_small(:,1)==B & CPT_flow_small(:,3)<0);
        % Consequently, there is a flux from B to i, if these compartments
        % are in the same order within a row and the flux in the third
        % column is positive OR the order of B and i are switched and the
        % flux in the third column is negative. This condition is met in
        % one logical array array_B_i
        array_B_i = (CPT_flow_small(:,1)==B & CPT_flow_small(:,2)==i & CPT_flow_small(:,3)>=0) | (CPT_flow_small(:,2)==B & CPT_flow_small(:,1)==i & CPT_flow_small(:,3)<0);
        
        % calculating the flux from i to B and B to i with the absolute
        % value since the order in CPT1_CPT2_flow contains the information
        % of the flow direction
        flow_i_B = sum(abs(CPT_flow_small(array_i_B,3)));
        flow_B_i = sum(abs(CPT_flow_small(array_B_i,3)));
        
        % saving the information within CPT1_CPT2_flow
        CPT1_CPT2_flow(row,1) = i;
        CPT1_CPT2_flow(row,2) = B;
        CPT1_CPT2_flow(row,3) = flow_i_B;
        CPT1_CPT2_flow(row+1,1) = B;
        CPT1_CPT2_flow(row+1,2) = i;
        CPT1_CPT2_flow(row+1,3) = flow_B_i;
        
        % increase the row counter
        row = row + 2;
        overall_time = toc;
        delta_t = overall_time - start;
        % writing the information into the fileID
        fprintf(fileID,'Overall time [s]: %d\n',overall_time);
        fprintf(fileID,'Time for calculating the exchange stream [s]: %d\n',delta_t);
        fprintf(fileID,'Mass flow from CPT %d to CPT %d [kg/s]: %d\n',i,B,flow_i_B);
        fprintf(fileID,'Mass flow from CPT %d to CPT %d [kg/s]: %d\n',B,i,flow_B_i);       
    end
    % delete all used rows in CPT_flow_copy in order to avoid doubled calculation
    CPT_flow_copy(sum(CPT_flow_needed,2)>0,:)=[];
end

% same caluclation method for the inlet and the outlet

% determine the compartments, which are connected to the inlet
inlet_unique = unique(CPT1_inlet_flow(:,1));
% determine the number of compartments, which are connected to the inlet
n_inlet_CPT = size(inlet_unique,1);

CPT1_inlet_2 = zeros(n_inlet_CPT*2,2);
row = 1;
for i=1:n_inlet_CPT
    start = toc;
    B = inlet_unique(i);
    array_inlet_B = (CPT1_inlet_flow(:,1)==-1 & CPT1_inlet_flow(:,2)==B & CPT1_inlet_flow(:,3)>=0) | (CPT1_inlet_flow(:,1)==B & CPT1_inlet_flow(:,2)==-1 & CPT1_inlet_flow(:,3)<0);
    array_B_inlet = (CPT1_inlet_flow(:,1)==B & CPT1_inlet_flow(:,2)==-1 & CPT1_inlet_flow(:,3)>=0) | (CPT1_inlet_flow(:,1)==-1 & CPT1_inlet_flow(:,2)==B & CPT1_inlet_flow(:,3)<0);
    flow_inlet_B = sum(abs(CPT1_inlet_flow(array_inlet_B,3)));
    flow_B_inlet = sum(abs(CPT1_inlet_flow(array_B_inlet,3)));
    CPT1_inlet_2(row,1) = -1;
    CPT1_inlet_2(row,2) = B;
    CPT1_inlet_2(row,3) = flow_inlet_B;
    CPT1_inlet_2(row+1,1) = B;
    CPT1_inlet_2(row+1,2) = -1;
    CPT1_inlet_2(row+1,3) = flow_B_inlet;
    row = row + 2;
    overall_time = toc;
    delta_t = overall_time - start;
    
    fprintf(fileID,'Overall time [s]: %d\n',overall_time);
    fprintf(fileID,'Time for calculating the exchange stream [s]: %d\n',delta_t);
    fprintf(fileID,'Mass flow from inlet to CPT %d [kg/s]: %d\n',B,flow_inlet_B);
    fprintf(fileID,'Mass flow von CPT %d to inlet [kg/s]: %d\n',B,flow_B_inlet);
 end

CPT1_inlet_2(CPT1_inlet_2(:,1) == 0,:) = [];

%---------------------------

% determine the compartments, which are connected to the outlet
outlet_unique = unique(CPT1_outlet(:,1));
% determine the number of compartments, which are connected to the outlet
n_outlet_CPT = size(outlet_unique,1);

CPT1_outlet_2 = zeros(n_outlet_CPT*2,2);
row = 1;
for i=1:n_outlet_CPT
    start = toc;
    B = outlet_unique(i);
    array_outlet_B = (CPT1_outlet(:,1)==-2 & CPT1_outlet(:,2)==B & CPT1_outlet(:,3)>=0) | (CPT1_outlet(:,1)==B & CPT1_outlet(:,2)==-2 & CPT1_outlet(:,3)<0);
    array_B_outlet = (CPT1_outlet(:,1)==B & CPT1_outlet(:,2)==-2 & CPT1_outlet(:,3)>=0) | (CPT1_outlet(:,1)==-2 & CPT1_outlet(:,2)==B & CPT1_outlet(:,3)<0);
    flow_outlet_B = sum(abs(CPT1_outlet(array_outlet_B,3)));
    flow_B_outlet = sum(abs(CPT1_outlet(array_B_outlet,3)));
    CPT1_outlet_2(row,1) = -2;
    CPT1_outlet_2(row,2) = B;
    CPT1_outlet_2(row,3) = flow_outlet_B;
    CPT1_outlet_2(row+1,1) = B;
    CPT1_outlet_2(row+1,2) = -2;
    CPT1_outlet_2(row+1,3) = flow_B_outlet;
    row = row + 2;
    overall_time = toc;
    delta_t = overall_time - start;
    
    fprintf(fileID,'Overall time [s]: %d\n',overall_time);
    fprintf(fileID,'Time for calculating the exchange stream [s]: %d\n',delta_t);
    fprintf(fileID,'Mass flow from outlet to CPT %d [kg/s]: %d\n',B,flow_outlet_B);
    fprintf(fileID,'Mass flow von CPT %d to outlet [kg/s]: %d\n',B,flow_B_outlet);
 end

CPT1_outlet_2(CPT1_outlet_2(:,1) == 0,:) = [];

%-----------------------------

% delete all rows with 0 within CPT1_CPT2_flow in the first column since these entries were not
% overwritten
CPT1_CPT2_flow(CPT1_CPT2_flow(:,1) == 0,:) = [];

% create one matrix with all relevant flows
CPT1_CPT2_flow = [CPT1_CPT2_flow;CPT1_inlet_2;CPT1_outlet_2];

start = toc;

% create a matrix to check for continuity of every compartment, the inlet
% and the outlet
balance = zeros(n_CPT+2,4);

% start at -2 in order to consider both the inlet and the outlet
for i = -2:n_CPT
    % if the compartment is in the first column, the flow in the third
    % column leaves this compartmentm. If the compartment is in the second
    % column, the flow enters the compartment
    outflow = sum(CPT1_CPT2_flow(CPT1_CPT2_flow(:,1)==i,3));
    inflow = sum(CPT1_CPT2_flow(CPT1_CPT2_flow(:,2)==i,3));
    balance(i+3,1) = i;
    balance(i+3,2) = inflow;
    balance(i+3,3) = outflow;
    balance(i+3,4) = inflow - outflow;
end

% delete all rows with a 0 within balance since there is no compartment 0
balance(balance(:,1)==0,:) = [];

overall_time = toc;
delta_t = overall_time - start;

fprintf(fileID,'Time for balance calculation [s]: %d\n',delta_t);

start = toc;

% delete all rows within CPT1_CPT2_flow with zero flow in the third column
% to reduce the number of entries.
CPT1_CPT2_flow(CPT1_CPT2_flow(:,3)==0,:)=[];

%%
% creating a Flow_Matrix, which contains all exchange streams between all
% compmartments and the inlet and outlet. A row represenets all in- and
% outcoming streams for one compartment. All streams, which leave one
% compartment x, are listed within column x. The inlet and outlet are
% treated as a compartment as well within the matrix and are represented by
% the index N_CPT + 1 and n_CPT + 2 respectively.

% since the Flow_matrix is a square matrix, its memory space increases by
% the square of number of compartments. In order to reduce the needed
% memory space, a sparese matrix is created, which only contains the
% indexes of non-zero elements within this matrix and its value. The
% Flow_Matrix is built like the following:
%       1   2   3 ... -1 -2 (#CPT, Inlet, Outlet)  
%  1
%  2
%  3
%  .
%  .
%  .
%  -1
%  -2

% For example: entry (2,3) is the flux from compartment 3 into compartment
% 2. Accordingly, entry (3,2) is the flux from compartment 2 to compartment
% 3. All entries within a row x, except the element on the diagonal, are
% the incoming fluxes for compartment x. All entries within a column y,
% except the element on the diagonal, are the outgoing fluxes for
% compartnemt y (Flow_Matrix(to, from)).

% since CPT1_CPT2_flow already is a list of all fluxes (column 3) from a
% compartment (column 1) to a comparmtent (column 2), the matrix is
% transformed to a sparse matrix. Before that, the entries for the inlet
% and outlet (-1 -2) must be changed to positive values for the index.
% Since the inlet is in the forelast row/column and the outlet is in the
% last row/column, their indexes can be replaced by n_CPT + 1 and n_CPT + 2
% respectively;

CPT1_CPT2_flow(CPT1_CPT2_flow(:,1:2) == -1) = n_CPT + 1;
CPT1_CPT2_flow(CPT1_CPT2_flow(:,1:2) == -2) = n_CPT + 2;

% creating the matrix with the sparse function. The last two input
% arguments defines the matrix as a square matrix with the number of
% compartments + the inlet and outlet.

% sparse(to, from, flow, n_CPT + 2, n_CPT + 2)
Flow_Matrix = sparse(CPT1_CPT2_flow(:,2),CPT1_CPT2_flow(:,1),CPT1_CPT2_flow(:,3),n_CPT+2,n_CPT+2);

% the matrix does not contain any entries on the diagoanl yet, which
% represent the total ouflow of every compartment. Due to the mass balance
% for every compartment and the assumption that the density is constant,
% the overall outgoing flux for a compartment is equal to the sum of all
% ingoing fluxes. This is represented by the sum along every row. Since the
% flow is an outgoing flux, the sign is negative

D = diag(-sum(Flow_Matrix,2));
Flow_Matrix = Flow_Matrix + D;

% alle fluxes are mass flow rates. Because the balance for a CSTR is based
% on a volumetric change, the fluxes must be changed to volume flow rates
% (./rho) and every row must be "scaled" with the corresponding compartment
% volume (./Volume_all). The inlet and outlet are associated with a volume
% of 1. Thus, only the last two rows in the matrix contains the exact
% volume flow rates.

Volume_all = [CPT_V(:,2);1;1];
Flow_Normalized = Flow_Matrix./Volume_all./rho;

% Flows into the outlet
Flow_Outlet = Flow_Normalized(end,1:end-2);
% Volume of the domain
V_all = sum(CPT_V(:,2));
% Volume flow through the domain/outlet
Volume_flow = sum(Flow_Outlet);
% hydrodynamic residence time
tau = V_all/Volume_flow;

% creation of the simplified matrix for the TIS model
Flow_Matrix_TIS = -eye(n_SLICES+2, n_SLICES+2)*Volume_flow; % n_SLICES + 2 for the inlet and outlet
% writing the elements below the diagonal
Flow_Matrix_TIS(2:(n_SLICES+2+1):(n_SLICES+2)^2) = Volume_flow;
% manipulate entries for the inlet
Flow_Matrix_TIS(end,end-1) = 0;
Flow_Matrix_TIS(1,end-1) = Volume_flow;
Flow_Matrix_TIS(end-1,end-2) = 0;
Flow_Matrix_TIS(end,end-2) = Volume_flow;

Flow_Outlet_TIS = Flow_Matrix_TIS(end,1:end-2);

V_CSTR = sum(CPT_V(:,2))/n_SLICES;
tau_CSTR = Volume_flow/V_CSTR;
Flow_Matrix_TIS_Normalized = Flow_Matrix_TIS/V_CSTR;

overall_time = toc;
delta_t = overall_time - start;

fprintf(fileID,'Time for creating the flow matrix [s]: %d\n',delta_t);
fprintf(fileID,'Time for calculating the streams [s]: %d\n',overall_time);
fclose(fileID);

save(save_string,'c0_M1','c0_c1_flow','c0_inlet_flow','c0_outlet_flow',...
    'CPT_V','c0_x_y_z_V','Barycenter','Vol_Avg_M1','CPT','delta_x',...
    'n_SLICES','n_LOCAL','CPT_slice','Flow_Normalized','CPT_V','rho','balance',...
    'Flow_Outlet','V_all','Volume_flow','tau','n_CELLS','n_CPT','Flow_Matrix_TIS',...
    'Flow_Matrix_TIS_Normalized','Flow_Outlet_TIS','V_CSTR','tau_CSTR');    

end

