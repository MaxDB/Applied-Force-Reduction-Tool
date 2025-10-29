function Static_Data = one_mode_rom(system_name,energy_limit,initial_modes,Static_Opts)
Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
Static_Data = Static_Dataset(Model);
Static_Data.save_data;
end