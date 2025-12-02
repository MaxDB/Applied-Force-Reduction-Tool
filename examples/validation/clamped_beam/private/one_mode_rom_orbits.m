function Dyn_Data = one_mode_rom_orbits(system_name,state)
Dyn_Data = initalise_dynamic_data(system_name);

%-------------------------------------------------------------------------%
%--------- Continuation Settings ---------%
Continuation_Opts.initial_inc = 5e-2;
Continuation_Opts.max_inc = 5e-2;
Continuation_Opts.min_inc = 5e-3;
Continuation_Opts.forward_steps = 200;
Continuation_Opts.backward_steps = 0;
Continuation_Opts.collocation_degree = 6;
Continuation_Opts.initial_discretisation_num = 20;
Continuation_Opts.min_discretisation_num = 20;
%-----------------------------------------%

switch state
    case 1
        Dyn_Data = Dyn_Data.add_backbone(1,"opts",Continuation_Opts);
    case 2
        Continuation_Opts.initial_inc = 2e-3;
        Continuation_Opts.min_inc = 1e-3;
        Continuation_Opts.max_inc = 2e-3;
        Continuation_Opts.forward_steps = 0;
        Continuation_Opts.backward_steps = 250;
        %-----------------------------------------
        Dyn_Data = Dyn_Data.add_orbits(1,[7,11],"opts",Continuation_Opts);
end
end