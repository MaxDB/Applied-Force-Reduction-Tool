function Static_Data = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts)
Verification_Opts.verification_algorithm = "sep_to_edge";

Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model,Verification_Opts);
Static_Data.save_data;
end