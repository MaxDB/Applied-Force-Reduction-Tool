clear
close all

num_iterations = 10;

%--------- Software Settings ---------%
set_logging_level(3)
set_visualisation_level(0)
%-------------------------------------%

%--------- System Settings ---------%
system_name = "clamped_beam";
energy_limit = 0.01;
initial_modes = [1];
%-----------------------------------%

%--------- Static Solver Settings ---------%
Static_Opts.max_parallel_jobs = 8; %be careful!
%------------------------------------------%
if isempty(gcp('nocreate'))
    parpool;
end

rom_one_base = zeros(1,num_iterations);
rom_one_validation = zeros(1,num_iterations);

for iCount = 1:num_iterations
    delete_static_data(system_name+"_"+initial_modes);
    base_time_start = tic;
    Static_Opts.additional_data = "none";
    Static_Data_Base = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts);
    rom_one_base(iCount) = toc(base_time_start);
    
    delete_static_data(system_name+"_"+initial_modes);
    validation_time_start = tic;
    Static_Opts.additional_data = "stiffness";
    Static_Data_Validaition = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts);
    rom_one_validation(iCount) = toc(validation_time_start);
end




%-----------------------------------
function Static_Data = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts)
Verification_Opts.verification_algorithm = "sep_to_edge";

delete_cache(system_name,"force",energy_limit)
delete_cache(system_name,"matrices")

Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end