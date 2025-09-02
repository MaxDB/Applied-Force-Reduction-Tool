function Model = create_model(system_name,energy_limit,initial_modes,Static_Opts)
Model = Dynamic_System(system_name,energy_limit,initial_modes,"static_opts",Static_Opts);
end