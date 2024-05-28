clear
compare_solutions("h_oscillator_1","h_oscillator_12",1,1,"energy",1)
compare_solutions("h_oscillator_1","h_oscillator_12",1,1,"amplitude",1)

compare_solutions("h_beam_1","h_beam_13",1,1,"amplitude",0)
compare_solutions("h_beam_1","h_beam_13",2,1,"amplitude",1)

plot_h_predicition("exhaust_1","energy",1);
plot_h_predicition("exhaust_1","amplitude",1);
compare_solutions("exhaust_17","exhaust_157",2,1,"energy",0)

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