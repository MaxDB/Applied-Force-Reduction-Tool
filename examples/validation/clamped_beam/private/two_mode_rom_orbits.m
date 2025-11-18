function Dyn_Data = two_mode_rom_orbits(system_name)
Continuation_Opts.initial_inc = 1e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.min_inc = 1e-3;
Continuation_Opts.forward_steps = 500;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.max_discretisation_num = 250;
Continuation_Opts.min_discretisation_num = 20;
Continuation_Opts.collocation_degree = 6;
%-----------------------------------------%
Dyn_Data = initalise_dynamic_data(system_name);

Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
end