function [] = FUNC_CLUSTERING(n_SLICES,n_LOCAL)
% n_SLICES: number of slices for the domain
% n_LOCAL: divisions of local Mean-Age within a slice

x_max = 55.6e-3;    % x-position of the outlet in m
x_min = -2e-3;      % x-position of the inlet in m
delta_x = (x_max-x_min)/n_SLICES;   % width of slice in m

% preparing the variable strings for loading, saving, the diary and the
% interpolation file
load_string = 'DATA_PREPARED_RE_1.mat';
save_string = strcat('CLUSTERED_DATA_RE_1_SLICES_',num2str(n_SLICES),'_DELTA_',num2str(n_LOCAL),'.mat');
diary_string = strcat('Diary_slices_RE_1_',num2str(n_SLICES),'_delta_',num2str(n_LOCAL),'_clustering.txt');
interpolation_string = strcat('INTERPOLATION_RE_1_SLICES_',num2str(n_SLICES),'_DELTA_',num2str(n_LOCAL),'_CPT_NUMBER.txt');

% load the relevant data for the clusteirng process
load(load_string,'c0_M1','c0_c1_flow','c0_inlet_flow','c0_outlet_flow','c0_x_y_z_V');


% create a fileID for writing the diary
fileID = fopen(diary_string,'w');
fprintf(fileID,('\n-----START: Agglomeration-----\n'));

% set the counter for the compartments to 0 for the beginning
i=0;

% create a vector CPT, with the the dimension of n_cells and NaN as entries
% as a dummy. Entry i equals to cell i-1 fromm FLUENT. CPT will contain the
% compartment number for every cell
n_CELLS=size(c0_M1,1);
CPT=zeros(n_CELLS,1)*NaN;

% create a matrix CPT_slice with the CPT number in the first column and the
% slice number in second column in order to save whicht compartment belongs
% to which slice.
CPT_slice = zeros(n_CELLS,2);

% create a matrix with the slice number in the first column, the volume
% average, variance and CoV within this slice w.r.t to M1 in order to
% compare the CFD and the CPT model regarding the radial resolution
Slice_M1_var_CoV = zeros(n_SLICES,4);

% create a matrix with the needed entries for the key numbers
M1_x_V = [c0_M1(:,2) c0_x_y_z_V(:,[2 5])];

% y is the number of already used cells
n_cells_overall=0;

% start of the time measurement
t=0;

% create a matrix in order to save the total volume of every compartnemt
% (#comp V_comp). Since it is nor possibly to know the number of
% compartments a priori, the matrix is created as a zero-matrix with
% n_cells as rows. Entries with a 0 will be deleted afterwards.
CPT_V = zeros(n_CELLS,2);

% create a matrix in order to save the barycenter of every compartment (x_s
% y_s z_s)
Barycenter = zeros(n_CELLS,3);

% create a mtrix in order to save the volume averaged value of the
% clustered cells for every compartment (#comp Vavg_M1)
Vol_Avg_M1 = zeros(n_CELLS,2);

% cut the domain into slices and call FUNC_CLUSTERING for every slice
for j=1:n_SLICES
    % max and min x-position for the slice
    x1 = x_min + (j-1)*delta_x;
    x2 = xmin + j*delta_x;
    % possible_cells contains only cells within the specified x-range
    possible_cells = c0_x_y_z_V(c0_x_y_z_V(:,2)>x1 & c0_x_y_z_V(:,2)<= x2,1);
    % c0_M1_small is in the first instance a locigal array, whicht contains
    % the value true, if the cells in c0_M1 are also in possible_cells.
    % Afterwards, these cells are extracted from c0_M1 and c0_M1_small is a
    % matrix
    c0_M1_small = ismember(c0_M1(:,1),possible_cells);
    c0_M1_small = c0_M1(c0_M1_small,:);
    % c0_c1_flow_small is in the first instance a logical array, which
    % contains the value true, if the cells in c0_c1_flow are also in
    % possible_cells. Afterwards, these cells are extracted from c0_c1_flow
    % and c0_c1_flow_small is a matrix
    c0_c1_flow_small = ismember(c0_c1_flow(:,1),possible_cells) & ismember(c0_c1_flow(:,2),possible_cells);
    c0_c1_flow_small = c0_c1_flow(c0_c1_flow_small,:);
    % defining the delta for the algorithm in the slice
    delta = (max(c0_M1_small(:,2))-min(c0_M1_small(:,2)))/n_LOCAL;
    % c0_M1_small will be overwritten in the second column with a NaN, if
    % the cell has been used for a compartment. As long as there is a cell,
    % its M1 has not been overwritten yet, the logical statement will
    % contain at least one true entry so that the while loop will continue
    while find(isnan(c0_M1_small(:,2))~=true)>0

    % starting the time measurement
    tic

    % increase the compartment counter
    i=i+1;
    
    % calling the function for the clustering process of one compartment
    % Input: - c0_M1_small: all cells with its parameter within the slice -
    % c0_c1_flow_small: all connectivites and flows within the slice -
    % delta: specifed deviation for the starting cell Output: - C: a vector
    % which contains all cell indices for the generated compartment
    % (adjusted to the MATLAB enumeration) - M1_max: the value of M1 of the
    % starting cell
    [C,M1_max]=FUNC_CLUSTERN_FIXED_DELTA_SLICES(c0_M1_small, c0_c1_flow_small,delta);

    % CPT is overwritten with the global compartment number i at every cell
    % ebtry
    CPT(C)=i;
    
    % create a second C vector with 1 at every position for the current
    % compartment number and 0 for all other entries for a fast calculation
    % of the compartment volume, barycenter etc.
    C2 = zeros(n_CELLS,1);
    C2(C) = 1;

    % calculation of the compartment related properties
    volumen = c0_x_y_z_V(:,5);
    CPT_V(i,1)=i;
    CPT_V(i,2)=sum(volumen.*C2,1);
    Barycenter(i,1) = sum(c0_x_y_z_V(:,2).*c0_x_y_z_V (:,5).*C2)/CPT_V(i,2);
    Barycenter(i,2) = sum(c0_x_y_z_V(:,3).*c0_x_y_z_V (:,5).*C2)/CPT_V(i,2);
    Barycenter(i,3) = sum(c0_x_y_z_V(:,4).*c0_x_y_z_V (:,5).*C2)/CPT_V(i,2);
    Vol_Avg_M1(i,1) = i;
    Vol_Avg_M1(i,2) = sum(c0_M1(C2==1,2).*c0_x_y_z_V(C2==1,5))/sum(c0_x_y_z_V(C2==1,5));
    CPT_slice(i,1) = i;
    CPT_slice(i,2) = j;
    % all cells, which were used for the current compartment, are
    % overwritten in c0_M1_small with NaN so that they are not considered
    % for the next iteration
    c0_M1_small(ismember(c0_M1_small(:,1),C-1),2)=NaN;

    % end of time measurement for the current compartment
    z=toc;
    t=t+z;

    % number of used cells for the current compartment
    n_cells_cpt = sum(C2);

    % number of overall used cells
    n_cells_overall = n_cells_overall+sum(C2);

    % number of free cells
    n_cells_left = n_CELLS - n_cells_overall;

    % writing the infomration into fileID
    fprintf(fileID,'Compartment: %d \n',i);

    fprintf(fileID,'Value of seeding cell [s]: %d \n',M1_max);
    fprintf(fileID,'Used deviation [s]: %d \n',delta);
    fprintf(fileID,'Number of cells in compartment: %d \n',n_cells_cpt);
    fprintf(fileID,'Volume of compartment [m^3]: %d \n',CPT_V(i,2));
    fprintf(fileID,'Volume averaged Mean-Age [s]: %d \n',Vol_Avg_M1(i,2));
    fprintf(fileID,'Number of overall used cells: %d \n',n_cells_overall);
    fprintf(fileID,'Number of free cells: %d \n',n_cells_left);
    fprintf(fileID,'Overall time [s]: %d\n',t);
    fprintf(fileID,'Time for compartment generation [s]: %d\n ',z);
    end
    % calculate the key numbers for the slice
    M1_V_sub = M1_x_V(M1_x_V(:,2)>x1 & M1_x_V(:,2)<=x2,[1 3]);
    M1_Vavg = sum(M1_V_sub(:,1).*M1_V_sub(:,2))/sum(M1_V_sub(:,2));
    var = sum((M1_Vavg - M1_V_sub(:,1)).^2.*M1_V_sub(:,2))/sum(M1_V_sub(:,2));
    CoV = sqrt(var)/M1_Vavg;
    Slice_M1_var_CoV(j,1) = j;
    Slice_M1_var_CoV(j,2) = M1_Vavg;
    Slice_M1_var_CoV(j,3) = var;
    Slice_M1_var_CoV(j,4) = CoV;
end
fprintf(fileID,'Compartmentbildung beendet.\n');
fclose(fileID);

% determine the mean, variance and CoV of M1 within every slice for a
% comparision with the compartmentmodel

% delete all unnecessary entries
CPT_V(CPT_V(:,1)==0,:) = [];
Barycenter(Barycenter(:,1)==0,:) = [];
Vol_Avg_M1(Vol_Avg_M1(:,1)==0,:) = [];
CPT_slice(CPT_slice(:,1)==0,:) = [];

% writing interpolation file for FLUENT to show different compartments.
% Export the commpartment number as an additional UDS. Inputmatrix must
% have the form : x,y,z,#CPT
Inputmatrix = zeros(size(CPT,1), 3 + 1);
Inputmatrix(:,1:3) = c0_x_y_z_V(:,2:4);
n_CPT = max(CPT);
for i = 1:n_CPT
    Inputmatrix(CPT==i,end) = CPT(CPT==i);
end

exporttofluent(Inputmatrix,{'uds-6'},interpolation_string);

save(save_string,'c0_M1','c0_c1_flow','c0_inlet_flow','c0_outlet_flow',...
    'CPT_V','c0_x_y_z_V','Barycenter','Vol_Avg_M1','CPT','delta_x',...
    'n_SLICES','n_LOCAL','n_CPT','n_CELLS','CPT_slice','Slice_M1_var_CoV')
end