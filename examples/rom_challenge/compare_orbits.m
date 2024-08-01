clear
close all
system_name = "exhaust_1567";
Dyn_Data = initalise_dynamic_data(system_name);

num_sols = size(Dyn_Data,1);
orbit_index = zeros(0,2);
for iSol = 1:num_sols
    sol_orbit_index = Dyn_Data.get_special_point(iSol,"X");
    num_orbits = size(sol_orbit_index,1);
    sol_id = iSol*ones(num_orbits,1);
    orbit_index = [orbit_index;[sol_id,sol_orbit_index]]; %#ok<AGROW>
end


% orbit_groups = {[1:4,11:13];5:7;8:10};
orbit_groups = {1:4,[5,10,11],[6,9,12],[7,8]};
num_groups = length(orbit_groups);
colour_num = [4,3,2];
colour_groups = {colour_num(1)*ones(1,4),colour_num([3,1,2]),colour_num([3,1,2]),colour_num([3,1])};

for iGroup = 1:num_groups
    ax = [];
    plot_orbit_index = orbit_index(orbit_groups{iGroup},:);

    num_orbits = size(plot_orbit_index,1);
    for iOrbit = 1:num_orbits
        orbit_colour = colour_groups{iGroup}(iOrbit);
        ax = plot_orbit(Dyn_Data,"physical displacement",plot_orbit_index(iOrbit,1),plot_orbit_index(iOrbit,2),"axes",ax,"colour",orbit_colour);
    end
end