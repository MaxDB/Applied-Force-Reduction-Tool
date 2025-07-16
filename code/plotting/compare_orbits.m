function ax = compare_orbits(type,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

ax = [];
plot_legend = 1;
stability = 1;
colour_num = [];
tag = "";

dyn_data_names = cell(0,1);
orbit_indices = cell(0,1);
system_counter = 0;
for arg_counter = 1:num_args/2
    if class(keyword_args{arg_counter}) == "Dynamic_Dataset"
        system_counter = system_counter + 1;
        dyn_data_names{1,system_counter} = keyword_args{arg_counter};
        orbit_indices{1,system_counter} = keyword_values{arg_counter};
        continue
    end
    switch keyword_args{arg_counter}
        case "validation"
            validation = keyword_values{arg_counter};
        case "axes"
            ax = keyword_values{arg_counter};
        case "legend"
            plot_legend = keyword_values{arg_counter};
        case {"color","colour"}
            colour_num = keyword_values{arg_counter};
        case {"stability"}
            stability = keyword_values{arg_counter};
        case "tag"
            tag = keyword_values{arg_counter}; 
        otherwise
            system_counter = system_counter + 1;
            dyn_data_names{1,system_counter} = keyword_args{arg_counter};
            orbit_indices{1,system_counter} = keyword_values{arg_counter};  
    end
end
%-----------------------------


for iSystem = 1:system_counter
    system = dyn_data_names{iSystem};
    Dyn_Data = initalise_dynamic_data(system);
    
    orbit_index = orbit_indices{iSystem};
    switch class(orbit_index)
        case "double"
            solution_num = orbit_index(:,1);
            orbit_num = orbit_index(:,2);
        case "cell"
            num_solutions = size(orbit_index,1);
            solution_num = zeros(0,1);
            orbit_num = zeros(0,1);
            orbit_counter = 0;
            for iSol = 1:num_solutions
                sol_num =  orbit_index{iSol,1};
                orbits = orbit_index{iSol,2};
                if isstring(orbits) && orbits == "all"
                    Sol = Dyn_Data.load_solution(sol_num);
                    orbits = 1:Sol.num_orbits;
                end
                num_orbits = length(orbits);
                for iOrbit = 1:num_orbits
                    orbit_counter = orbit_counter + 1;
                    
                    solution_num(orbit_counter,1) = sol_num; %#ok<*AGROW>
                    orbit_num(orbit_counter,1) = orbits(iOrbit);
                end
                
            end
    end

    if isempty(colour_num)
        orbit_colour = iSystem;
    else
        orbit_colour = colour_num;
        if length(orbit_colour) > 1
            orbit_colour = orbit_colour(iSystem);
        end
    end
    
    num_orbits = length(orbit_num);
    for iOrbit = 1:num_orbits
        ax = plot_orbit(Dyn_Data,type,solution_num(iOrbit),orbit_num(iOrbit),"axes",ax,"colour",orbit_colour,"stability",stability,"tag",tag);
    end
end

end
