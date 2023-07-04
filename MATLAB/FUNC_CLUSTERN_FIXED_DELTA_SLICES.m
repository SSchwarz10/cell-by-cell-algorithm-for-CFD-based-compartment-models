function [C,p_max]=FUNC_CLUSTERN_FIXED_DELTA_SLICES(c0_P, c0_c1_flow, delta)
%% subfunction for clustering the CFD cells, which is called for every CPT within the function "FUNC_CLUSTERING"
% Input:
% - c0_P: matrix containing the cells and its parameter (c0,P)
% - c0_c1_flow_small: matrix with all connectivites and flows (c0,c1,P)
% - delta: specifed deviation for the starting cell
% Output:
% - C: vector which contains all cell indices for the generated
% compartment (adjusted to the MATLAB enumeration)
% - p_max: the value of M1 of the starting cell

% extract all cells of c0_P which are in the considered range
c0_P_small=c0_P((c0_P(:,2)>=(max(c0_P(:,2))-delta)),:);
% consider only the connectivity of cells
connectivity=c0_c1_flow(:,1:2);
% identification of the maximum and the position within c0_P_small. index
% IS NOT the cell number
[p_max,index]=max(c0_P_small(:,2));


% create a logical array, which contains true, if connectivity has cells in
% the firt AND second column, which are within c0_P_small. If a row
% contains only one cell of c0_P_small, this row is not considered since
% only one cell satisfy the condition
IS_MEMBER_connectivity  = ismember(connectivity(:,1),c0_P_small(:,1)) & ismember(connectivity(:,2),c0_P_small(:,1));

% Possible_Rows extract the relevant rows from connectivity. Since it is
% possible that the starting cell has no neighbored cell, which satisfy the
% condition, it is possible, that Possible_Rows is empty!
Possible_Rows = connectivity(IS_MEMBER_connectivity,:);

% after all possible cells and rows are identified, it must be checked,
% whether cells form a continuous cluster or not. Starting from the
% maximum, neighbored cells are clusters successively. If no more cells are
% added to the cluster, the process is stopped

% c_unique contains the starting cell (maximum)
c_unique = c0_P_small(index,1);

% since it is possible, that the starting cell has no neighbored cells, it
% is useful to check, whether there are neighbored cells or not only once.

% the procedure: a logical array is created, which identifies all rows in
% Possible_Rows, in which all cells within c_unique are in either the first
% column or the second column. If the starting cell is not part of the
% Possible_Rows matrix, the logical array only contains 'false/0' and
% c_unique_new is empty. If there are neighbored cells, c_unique_new
% contains both the starting cell and all neighbored cells, which satisfy
% the condition for the cluster. Within the while loop, c_unique is
% overwritten with c_unique_new and contains more and more cells with every
% iteration
logical_array = ismember(Possible_Rows(:,1),c_unique) | ismember(Possible_Rows(:,2),c_unique);
c_unique_new = unique(Possible_Rows(logical_array,:));

% check, whether c_unique_new is empty. If so, no clustering prcoess is
% needed
if numel(c_unique_new) ~= 0
    % If c_unique_new is not empty, the cluster process starts.
    stop = 1;
    while stop ~= 0
    logical_array = ismember(Possible_Rows(:,1),c_unique) | ismember(Possible_Rows(:,2),c_unique);
    c_unique_new = unique(Possible_Rows(logical_array,:));
    % c_unique_new contains all cells, which belong to one continuous
    % cluster, only once und grows with every iteration. If the number of
    % cells equals to the number of cells during the last iteration, the
    % clustering process is stopped. If not, c_unique is overwritten with
    % c_unique_new and a new iteration is performed
    if numel(c_unique) == numel(c_unique_new)
        stop = 0;
        break;
    else
        c_unique = c_unique_new;
    end
    end
end
% C contains all cells of the generated compartment with the enumeration of
% FLUENT. By adding 1, the index agrees with the enumeration of MATLAB
C = c_unique+1;
end