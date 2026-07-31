clear
close all
%---------------------------------
% 0.01    --> 7,000
% 0.005   --> 24,500
% 0.0044 is too big
% 0.003   --> 100,000
% 0.002   --> 300,000
% 0.00125 --> 1,000,000
% 0.001   --> 2,000,000

% dof ≈ 0.0458 * seed_size ^ -2.53
%-------------------------------
dof_range = [5e3,2e6];
num_points = 50;

dofs = logspace(log10(dof_range(1)),log10(dof_range(2)),num_points);
seed_sizes = (dofs/0.0488).^(-1/2.53);
%---
num_workers = 1;

%-----
set_logging_level(2)
set_visualisation_level(0)
system_name = "mems_arch";
energy_limit = 0.8;
modes = [1,6];

Calibration_Opts.calibration_scale_factor = 1;
%----
Static_Opts.max_parallel_jobs = num_workers;
Static_Opts.additional_data = "stiffness";
create_parallel_pool(num_workers);

%--
eig_values = zeros(2,num_points);
num_dof = zeros(1,num_points);

for iSeed = 1:num_points
    seed_size = seed_sizes(iSeed);
    
    delete_cache(system_name,"matrices")
    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);
    Model = Dynamic_System(system_name,0,[1,6]);
    eig_values(:,iSeed) = Model.reduced_eigenvalues;
end

frequency = sqrt(eig_values);
freq_size = 10.^floor(log10(frequency(:,end)));
freq_norm = frequency./freq_size;
freq_diff = freq_norm - freq_norm(:,end);


ax_style = @(ax) set_ax_style(ax);

figure
tiledlayout
ax = nexttile;
loglog(num_dof,freq_diff(1,:))
ax_style(ax)

ax = nexttile;
loglog(num_dof,freq_diff(2,:))
ax_style(ax)

Plot_Data.num_dof = num_dof;
Plot_Data.eig_values = eig_values;



function set_ax_style(ax)
    box(ax,"on")
    ylim([1e-4,1e-1])
end