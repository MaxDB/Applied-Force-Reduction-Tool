function [M_bc,K_bc,node_map] = matrices_abaqus(system_name)

setup_time_start = tic;
project_path = get_project_path;

G_ID =fopen("geometry\" + system_name+ "\" + system_name + ".inp");
geometry = textscan(G_ID,'%s','delimiter','\n');
fclose(G_ID);
geometry = geometry{1,1};

M_ID =fopen(project_path + "\fe_templates\abaqus\matrix_output.inp");
matrix_template = textscan(M_ID,'%s','delimiter','\n');
fclose(M_ID);
matrix_template = matrix_template{1,1};

input_file = [geometry;matrix_template];
file_name = "matrix";


%Write matrix extraction file
wid = fopen("temp\" + file_name + ".inp","w");
for k=1:numel(input_file)
    fprintf(wid,'%s\r\n',input_file{k,1});
end
fclose(wid);

setup_time = toc(setup_time_start);
log_message = sprintf("Matrix input file created: %.1f seconds" ,setup_time);
logger(log_message,3)

%Execute in cmd
abaqus_time_start = tic;

project_directory = pwd;
cd temp
[status,cmdout] = system("abaqus j=" + file_name); %#ok<ASGLU>

while ~isfile(file_name + ".dat")
    pause(0.1)
end
while isfile(file_name + ".lck")
    pause(0.1)
end
cd(project_directory)
abaqus_time = toc(abaqus_time_start);
log_message = sprintf("Abaqus matrix analysis complete: %.1f seconds" ,abaqus_time);
logger(log_message,3)


%Read mass and stiffness
data_processing_time_start = tic;

M_file = "temp\" + file_name + "_MASS1.mtx";
K_file = "temp\" + file_name + "_STIF1.mtx";

pre_M = load(M_file);
pre_K = load(K_file);

cd(project_directory)

% Apply boundary conditions
K = sparse(pre_K(:,1),pre_K(:,2),pre_K(:,3));
M = sparse(pre_M(:,1),pre_M(:,2),pre_M(:,3));

%%% Added for Xiao Xiao's beam
% zero stiffness but still appearing with sparse output?
% due to constraint?
zero_indicies = pre_K(pre_K(:,3) == 0,1:2);
ci_zero = zero_indicies(:,1);
%%


indicies = pre_K(pre_K(:,3) == 1e36,1:2);
ci = indicies(:,1);

K_bc = K;
K_bc([ci;ci_zero],:) = [];
K_bc(:,[ci;ci_zero]) = [];


M_bc = M;
M_bc([ci;ci_zero],:) = [];
M_bc(:,[ci;ci_zero]) = [];





% Caclulate node mapping
dof_bc = length(K_bc);
% dof = length(K) - length(ci_zero);
dof = length(K);

node_map = zeros(dof,2);
node_map(:,1) = (1:dof)';
node_map([ci;ci_zero],:) = [];
node_map(:,1) = node_map(:,1) - size(ci_zero,1);
node_map(:,2) = (1:dof_bc)';

%need to account for coupling

data_processing_time = toc(data_processing_time_start);
log_message = sprintf("Matrices processed: %.1f seconds" ,data_processing_time);
logger(log_message,3)


end