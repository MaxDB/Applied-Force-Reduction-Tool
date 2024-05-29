function [index_ranges,range_stability] = get_stability_boundaries(stability,bifurcations)
num_orbits = size(stability,2);
stability_transition = zeros(0,2);
for iOrbit = 1:(num_orbits-1)
    if stability(iOrbit) + stability(iOrbit+1) == 1 %stability transition
        stability_transition = [stability_transition;[iOrbit,iOrbit+1]]; %#ok<AGROW>
    end
end
bifurcation_types = fields(bifurcations);
num_bifurcation_types = size(bifurcation_types,1);
all_bifurcation_index = [];
for iType = 1:num_bifurcation_types
    bifurcation_type = bifurcation_types{iType,1};
    bifurcation_index = bifurcations.(bifurcation_type);
    if isempty(bifurcation_index)
        continue
    end
    all_bifurcation_index = [all_bifurcation_index;bifurcation_index]; %#ok<AGROW>
end

num_transitions = size(stability_transition,1);
num_ranges = num_transitions + 1;
index_ranges = zeros(num_ranges,2);
index_ranges(1,1) = 1;
index_ranges(end,2) = num_orbits;
for iTransition = 1:num_transitions
    transition = stability_transition(iTransition,:);
    is_bifurcation = ismember(transition,all_bifurcation_index);
    if all(is_bifurcation == 0)
        transition_orbit = transition(1);
    else
        transition_orbit = transition(is_bifurcation);
    end
    index_ranges(iTransition,2) = transition_orbit;
    index_ranges(iTransition+1,1) = transition_orbit;
end



range_stability = ones(1,num_ranges);
range_stability(2:2:num_ranges) = 0;
end