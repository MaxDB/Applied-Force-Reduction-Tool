clear
plot_backbone("exhaust_1567","physical amplitude",1);


plot_h_predicition("exhaust_1","energy",1);
plot_h_predicition("exhaust_1","amplitude",1);
compare_solutions("exhaust_157","exhaust_157",1,2,"energy",0)
compare_solutions("exhaust_1567","exhaust_1567",2,3,"energy",0)
compare_solutions("exhaust_1567","exhaust_1567",2,3,"physical amplitude",0)
plot_eigenvalues("exhaust_1567",2,411);

compare_solutions("exhaust_1","exhaust_1567",1,2,"energy",1)

compare_solutions("exhaust_17","exhaust_1716",1,1,"energy",0)

compare_solutions("exhaust_15716","exhaust_15716",3,4,"energy",0)
compare_solutions("exhaust_156716","exhaust_156716",1,3,"energy",0)

compare_solutions("energy","exhaust_156716",[1,3])
compare_solutions("energy","exhaust_15716",[1,3,4])

compare_solutions("energy","exhaust_1567","all")

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