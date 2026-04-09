function [displacement,energy,additional_data,additional_data_time] = read_abaqus_static_data(file_name,step_type,num_nodes,num_dimensions,step_incs)

num_dofs = num_nodes*num_dimensions;
num_steps = length(step_type);
num_static_steps = sum(step_type(:,1) == "static");

step_start_pattern = "STEP" + whitespacePattern + digitsPattern + whitespacePattern + "INCREMENT" + whitespacePattern + "1";
energy_pattern = "RECOVERABLE STRAIN ENERGY";
node_output_pattern = "N O D E   O U T P U T";
displacement_table_rows = join(["%u",repmat("%f",1,num_dimensions)]," ");
% table_pattern = digitsPattern + asManyOfPattern(whitespacePattern + digitsPattern,NUM_DIMENSIONS,NUM_DIMENSIONS); 

%------------------------------------------------------------------------%
static_incs = step_incs(step_type(:,1) == "static");
num_static_incs = sum(static_incs);
displacement = zeros(num_dofs,num_static_incs);
energy = zeros(1,num_static_incs);
%------------------------------------------------------------------------%
dat_ID = fopen("temp\" + file_name + ".dat");
abaqus_data = textscan(dat_ID,'%s','delimiter','\n');
fclose(dat_ID);
abaqus_data = abaqus_data{1,1};

%------------------------------------------------------------------------%
additional_data_time = 0;
if any(step_type == "perturbation")
    add_data_step_start = tic;
    loadcase_def = find(startsWith(abaqus_data,"*loadcase"),1);
    counter = 1;
    while startsWith(abaqus_data{loadcase_def+counter,1},"*loadcase")
        counter = counter + 1;
    end
    num_h_modes = counter;

    additional_data = zeros(num_dofs,num_h_modes,num_static_steps);
    additional_data_time = additional_data_time + toc(add_data_step_start);
else
    additional_data = [];
end

%------------------------------------------------------------------------%

step_start_lines = [find(matches(abaqus_data,step_start_pattern,'IgnoreCase',true));length(abaqus_data)];

static_step_counter = 0;
static_inc_counter = 0;
for iStep = 1:num_steps
    switch step_type(iStep,1)
        case "static"
            static_step_counter = static_step_counter + 1;
            step_span = step_start_lines(iStep):(step_start_lines(iStep+1)-1);
            step_data = abaqus_data(step_span,1);

            % Energy
            energy_line = step_data(startsWith(step_data,energy_pattern,'IgnoreCase',true));

            % Displacement
            output_start_line = find(startsWith(step_data,node_output_pattern,'IgnoreCase',true));
            num_incs = size(output_start_line,1);


            for iInc = 1:num_incs
                static_inc_counter = static_inc_counter + 1;
                %energy
                line_data = textscan(energy_line{iInc,1},"%*s %*s %*s %f");
                inc_energy = line_data{1,1};
                energy(1,static_inc_counter) = inc_energy;

                %displacement

                inc_start = output_start_line(iInc);
                if iInc < num_incs
                    inc_end = output_start_line(iInc+1)-1;
                else
                    inc_end = size(step_data,1);
                end
                node_data = step_data(inc_start:inc_end,1);


                step_displacement = zeros(num_nodes,num_dimensions);
                for iLine = 1:length(node_data)
                    line = node_data{iLine,1};
                    if isempty(line)
                        continue
                    end

                    line_data = textscan(line,displacement_table_rows);
                    if isempty(line_data{1,end})
                        continue
                    end

                    node_num = line_data{1,1};
                    step_displacement(node_num,:) = [line_data{1,2:end}];
                end
                %------------------%
                displacement(:,static_inc_counter) = reshape(step_displacement',num_dofs,1);
            end
       

        case "stiffness"
            add_data_step_start = tic;
            base_name = "temp\" + file_name + "_STIF";
            if num_static_steps == 0
                static_step_counter = static_step_counter + 1;
            end
            movefile(base_name + iStep + ".mtx", base_name + static_step_counter + ".mtx")
            additional_data_time = additional_data_time + toc(add_data_step_start);
        case "perturbation"
            add_data_step_start = tic;
            step_span = step_start_lines(iStep):(step_start_lines(iStep+1)-1);
            step_data = abaqus_data(step_span,1);
            if num_static_steps == 0
                static_step_counter = static_step_counter + 1;
            end

            % Displacement
            output_start_lines = find(startsWith(step_data,node_output_pattern,'IgnoreCase',true));
            output_start_lines(end+1) = length(step_data); %#ok<AGROW>
            for iMode = 1:num_h_modes
                output_start_line = output_start_lines(iMode);
                output_end_line = output_start_lines(iMode+1);
                node_data = step_data(output_start_line:output_end_line,1);

                loadcase_displacement = zeros(num_nodes,num_dimensions);
                for iLine = 1:length(node_data)
                    line = node_data{iLine,1};
                    if isempty(line)
                        continue
                    end

                    line_data = textscan(line,displacement_table_rows);
                    if isempty(line_data{1,end})
                        continue
                    end

                    node_num = line_data{1,1};
                    loadcase_displacement(node_num,:) = [line_data{1,2:end}];
                end
                % additional_data(:,iMode,static_step_counter) = reshape(loadcase_displacement',num_dofs,1) - displacement(:,static_step_counter);
                additional_data(:,iMode,static_step_counter) = reshape(loadcase_displacement',num_dofs,1);
            end
            % additional_data(:,:,static_step_counter) = additional_data(:,:,static_step_counter) - displacement(:,static_step_counter);
            additional_data_time = additional_data_time + toc(add_data_step_start);
            %------
    end

end

end