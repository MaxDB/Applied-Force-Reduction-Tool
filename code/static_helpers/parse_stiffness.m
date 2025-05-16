function Stiffness = parse_stiffness(step_list,file_name,dofs)
num_steps = length(step_list);
% num_steps = max(step_list);
base_name = "temp\" + file_name + "_STIF";
file_name_parts = split(file_name,"_");
if any(regexp(file_name_parts(end-1),'[0-9]'))
    job_id = double(file_name_parts(end-1));
else
    job_id = double(file_name_parts(end));
end

% dimensions = [dofs,dofs,num_steps];

Stiffness = Sparse_Stiffness_Pointer(dofs,"job",job_id);
Stiffness = Stiffness.add_loadcases(base_name,step_list);

% first_step = step_list(1);
% num_loadcases = size(step_list,2);
% 
% pre_K = load_mtx(base_name + first_step + ".mtx");
% 
% indicies = pre_K(pre_K(:,3) == 1e36,1:2);
% removed_section = indicies(:,1);
% num_rows = size(pre_K,1);
% remove_row = false(num_rows,1);
% for iRow = 1:num_rows
%     nz_index_row = pre_K(iRow,1:2);
%     remove_row(iRow) = any(ismember(nz_index_row,removed_section));
%     if remove_row(iRow)
%         continue
%     end
%     pre_K(iRow,1:2) = arrayfun(@(x) x - sum(gt(x,removed_section)),nz_index_row);
% end
% 
% 
% pre_K(remove_row,:) = [];
% 
% 
% Stiffness = Sparse_Stiffness_Array(dimensions,pre_K(:,1:2),1,num_loadcases);
% Stiffness{1} = pre_K;
% nonzero_data = Stiffness.nonzero_data;
% % 
% keep_row = ~remove_row;
% 
% if isempty(gcp('nocreate'))
%     for iStep = 2:num_steps
%         pre_K = load_mtx(base_name + iStep + ".mtx");
%         nonzero_data(:,iStep) = pre_K(keep_row,3);
%     end
% else
%     parfor iStep = 2:num_steps
%         pre_K = load_mtx(base_name + iStep + ".mtx");
%         nonzero_data(:,iStep) = pre_K(keep_row,3);
%         % pre_K = pre_K(:,3);
%         % pre_K(remove_row,:) = [];
%         % nonzero_data(:,iStep) = pre_K;
%     end
% end
% Stiffness.nonzero_data = nonzero_data;
end