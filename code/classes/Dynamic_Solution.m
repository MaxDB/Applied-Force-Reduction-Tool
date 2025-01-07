classdef Dynamic_Solution
    properties
        Continuation_Options
        Solution_Type

        orbit_labels
        num_orbits
    end
    
    methods
        function obj = Dynamic_Solution(Continuation_Opts)
            obj = obj.update_continuation_opts(Continuation_Opts);
        end
        %-----------------------------------------------------------------%
        function obj = update_continuation_opts(obj,Continuation_Opts)
            Default_Continuation_Opts = read_default_options("continuation");         
            obj.Continuation_Options = update_options(Default_Continuation_Opts,obj.Continuation_Options,Continuation_Opts);
        end 
        %-----------------------------------------------------------------%


        %-----------------------------------------------------------------%
        % Analysis 
        %-----------------------------------------------------------------%
        function Analysis_Output = analyse_solution(~,solution_num,Add_Output)
            solution_name = "temp\dynamic_sol_" + solution_num; 
            bd = coco_bd_read(solution_name);

            solution_data  = coco_bd_col(bd, {"po.period", "LAB", "TYPE", "eigs","ENERGY","MAX(x)","MIN(x)"});
            
            %------------------------%
            % Labels
            orb_labels = solution_data{1,2};
            
            %------------------------%
            % Frequency
            period = solution_data{1,1};
            orbit_frequency = 2*pi./period;
            
            %------------------------%
            % Energy
            orbit_energy = solution_data{1,5};
            
            %------------------------%
            % Bifurcations
            orbit_type = solution_data{1,3};
            orbit_bifurcation.bp_index = find(orbit_type == "BP");
            orbit_bifurcation.pd_index = find(orbit_type == "PD");
            orbit_bifurcation.ns_index = find(orbit_type == "TR");
            orbit_bifurcation.sn_index = find(orbit_type == "SN");
           
            %------------------------%
            % Stability
            eigs = solution_data{1,4};
            orbit_stability = max(abs(eigs),[],1);
            
            %------------------------%
            % Amplitude
             % sol = po_read_solution('',convertStringsToChars(solution_name),orbit_labels(iSol));
            num_modes = size(solution_data{1,6},1)/2;
            min_disp = solution_data{1,6}(1:num_modes,:);
            max_disp = solution_data{1,7}(1:num_modes,:);
            orbit_amplitude = abs(max_disp - min_disp)/2;

            switch Add_Output.output
                case "physical displacement"
                    additional_output  = coco_bd_col(bd, "DISP");
            end


            %------------------------%
            Analysis_Output.orbit_labels = orb_labels;
            Analysis_Output.num_orbits = size(orb_labels,2);
            Analysis_Output.frequency = orbit_frequency;
            Analysis_Output.energy = orbit_energy;
            Analysis_Output.amplitude = orbit_amplitude;
            Analysis_Output.stability = orbit_stability;
            Analysis_Output.bifurcations = orbit_bifurcation;

            switch Add_Output.output
                case "physical displacement"
                    Analysis_Output.additional_dynamic_output = additional_output;
            end
        end
        %-----------------------------------------------------------------%

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        function Solution = update_dynamic_solution(Solution,Analysis_Output)
            sol_properties = properties(Solution);
            output_fields = fieldnames(Analysis_Output);

            num_properties = size(sol_properties,1);
            for iProp = 1:num_properties
                sol_property = sol_properties{iProp,1};
                if ~ismember(sol_property,output_fields)
                    continue
                end
                Solution.(sol_property) = Analysis_Output.(sol_property);
            end
            
        end
    end

end