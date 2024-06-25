clear

compare_solutions("energy","exhaust_15717",[1,3])
compare_solutions("energy","exhaust_1567",[2,3])


compare_solutions("physical amplitude","exhaust_1",1,"exhaust_17",1,"exhaust_157",[2,3],"exhaust_1567",[2,3])

compare_solutions("physical amplitude","exhaust_1",1,"exhaust_1",4)

compare_solutions("physical amplitude","exhaust_1","all")

plot_h_predicition("exhaust_1567","energy",1);

compare_solutions("energy","beam_oscillator_1",1,"beam_oscillator_12",1)

% close all
% plot_backbone("h_beam_13","amplitude",1);
% plot_backbone("h_oscillator_1","amplitude",1);
% 
% 
% plot_orbit("h_oscillator_1","displacement",1,81);
% plot_orbit("h_oscillator_12","displacement",2,113);
% 
% plot_orbit("h_beam_13","displacement",1,237);
% 
% compare_solutions("h_oscillator_1","h_oscillator_12",1,1,"energy",1)

% plot_backbone("h_oscillator_12","amplitude",2);

% plot_h_predicition("h_beam_1","amplitude",1);