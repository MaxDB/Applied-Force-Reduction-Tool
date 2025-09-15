function [displacement,energy,additional_data,additional_data_time] = read_abaqus_fil_data(file_name,step_types,num_nodes,num_dimensions)

num_dofs = num_nodes*num_dimensions;
step_types(step_types == "zero") = [];
num_steps = length(step_types);
num_static_steps = sum(step_types(:,1) == "static");
num_stiffness_steps = sum(step_types(:,1) == "stiffness");

if num_stiffness_steps >0 && num_stiffness_steps   ~= num_static_steps
    error("Abaqus error")
end

num_output_steps = num_static_steps;
%------------------------------------------------------------------------%
displacement = zeros(num_dofs,num_static_steps);
energy = zeros(1,num_static_steps);
additional_data = [];
additional_data_time = 0;
%------------------------------------------------------------------------%

fil_id = fopen("temp\" + file_name + ".fil");
abaqus_data = textscan(fil_id,'%s','delimiter','*',"EndOfLine",""); %quicker to read all at once if there is enough memory
fclose(fil_id);



abaqus_data = regexprep(abaqus_data{1},'[\n\r]+',''); %update for nonwindows




Model_Data.num_dimensions = num_dimensions;
step_indicies = find(find_record_by_type(abaqus_data,["inc start","inc end"],Model_Data));

if diff(step_indicies((end-1):end)) == 1
    step_indicies(end) = [];
end
if diff(step_indicies(1:2)) == 1
    step_indicies(1) = [];
end


step_spans = zeros(num_output_steps,2);
step_spans(:,1) = step_indicies(1:2:(2*num_output_steps));
step_spans(:,2) = step_indicies(2:2:(2*num_output_steps));

output_counter = 0;
step_counter = 0;
for iStep = 1:num_steps
    step_type = step_types(iStep);
    switch step_type
        case "static"
            output_counter = output_counter + 1;
            step_counter = step_counter + 1;
            step_span = step_spans(output_counter,1):step_spans(output_counter,2);

            step_data = abaqus_data(step_span);
            % step_data = abaqus_data;

            energy_index = find_record_by_type(step_data,"energy",Model_Data);
            step_energy = read_energy(step_data(energy_index));

            disp_index = find_record_by_type(step_data,"displacement",Model_Data);

            step_displacement = read_displacements(step_data(disp_index),Model_Data);
            
            energy(output_counter ) = step_energy;
            displacement(:,output_counter )  = reshape(step_displacement',num_dofs,1);
        case "stiffness"
            step_counter = step_counter + 1;
            add_data_step_start = tic;
            base_name = "temp\" + file_name + "_STIF";
            if num_static_steps == 0
                static_step_counter = static_step_counter + 1;
            end
            file_found = 0;
            while ~file_found
                stiffness_name = base_name + step_counter + ".mtx";
                file_found = isfile(stiffness_name);
                if ~file_found
                    step_counter = step_counter + 1;
                end
            end
            movefile(stiffness_name, base_name + output_counter + ".mtx")
            additional_data_time = additional_data_time + toc(add_data_step_start);
    end
end


%------------------------------------------------------------------------%

end
%-------

function record_indicies = find_record_by_type(abaqus_data,record_type,Model_Data)
num_record_type = size(record_type,2);

for iRecord = 1:num_record_type
    record_size = [];
    switch record_type(iRecord)
        case "displacement" %displacement
            output_id = 101;
            num_dimensions = Model_Data.num_dimensions;
            record_size = 3 + num_dimensions; %num_entries, record_id, node, dof_1, ....
        case "energy" %increment start
            output_id = 1999;
        case "output request" %output request def
            output_id = 1911;
        case "inc start" %increment start and summary
            output_id = 2000;
        case "inc end" % increment end
            output_id = 2001;
        case "abaqus release"
            output_id = 1921;
        case "element def"
            output_id = 1900;
        case "node def" %node definition
            output_id = 1901;
        case "label ref"
            output_id = 1940;
        case "active dof"
            output_id = 1902;
            % case 1933 %??
            %     output_id = 1933;
            % case 1934 %??
            %     output_id = 1934;
            % case 1931 %??
            %     output_id = 1931;
            % case 1932 %??
            %     output_id = 1932;

            % case 1922 %??
            %     output_id = 1922;


            % case 1999 %?? energy?
            %     output_id = 1999;

        otherwise
            error("Unimplemented record type: '" + record_type + "'")

    end

    output_id = int2str(output_id);
    id_length = size(output_id,2);
    id_length = int2str(id_length);
    if size(id_length,2) == 1
        id_length = [' ',id_length];
    end


    output_pat = ['I',id_length,output_id];
    if isempty(record_size)
        start_pat = 'I' + wildcardPattern(2) + digitsPattern;
        output_pat = start_pat + output_pat;
    else
        record_size = int2str(record_size);
        size_length = int2str(size(record_size,2));
        if size(size_length,2) == 1
            size_length = [' ',size_length]; %#ok<AGROW>
        end
        start_pat = ['I',size_length,record_size];
        output_pat = [start_pat,output_pat]; %#ok<AGROW>
    end

    if iRecord == 1
        record_pat = output_pat;
    else
        record_pat = record_pat | output_pat;
    end
end

record_indicies = startsWith(abaqus_data,record_pat);

end
%-------
function displacement = read_displacements(disp_data,Model_Data)
num_dimensions = Model_Data.num_dimensions;
FLOAT_FIELD_WIDTH = 22;

disp_pattern = @(line_length) "%" + (line_length - (FLOAT_FIELD_WIDTH*num_dimensions + (num_dimensions -1))) + "*s " + ...
join(repmat("%" + FLOAT_FIELD_WIDTH + "f",[1,num_dimensions])," %*c ");

% tic
% [disp_table_1,test_2,test_3] = cellfun(@(line) scan_disp(line,disp_pattern),disp_data);
% toc
% % displacement = vertcat(disp_table{:});


num_lines = size(disp_data,1);
displacement = zeros(num_lines,num_dimensions);
for iLine = 1:num_lines
    data_line = disp_data{iLine};
    disp_line = textscan(data_line,disp_pattern(size(data_line,2)),1,'Delimiter','');
    displacement(iLine,:) = [disp_line{:}];
end


    % function [line_out_1,line_out_2,line_out_3]  = scan_disp(line_in,disp_pattern)
    %     line_out = textscan(line_in,disp_pattern(size(line_in,2)),1,'Delimiter','');
    %     % line_out = horzcat(line_out_cell{:});
    %     line_out_1 = line_out{1};
    %     line_out_2 = line_out{2};
    %     line_out_3 = line_out{3};
    % end

end
%-------
function energy = read_energy(energy_data)
line_length = size(energy_data{1},2);
line_data = textscan(energy_data{1},"%*" + (line_length - 16*23) +"s %*c %*22f %*c %22f",1,"Delimiter",'');
energy = line_data{1};
end

%-------
function record_entries = parse_data_record(data_record)
record_length = size(data_record,2);


record_index = 1;
entry_index = -1;
while record_index < record_length
    record_char = data_record(record_index);
    switch record_char
        case 'I' %integer: Innx...x - where there are nn digits: x...x 
            size_index = record_index + [1,2];
            num_digits = str2double(data_record(size_index));

            digits_index = record_index + 2 + (1:num_digits);
            int_value = str2double(data_record(digits_index));
            
            entry_index = entry_index + 1;
            if entry_index == 0
                num_entries = int_value - 1;
                record_entries = cell(num_entries,1);
            else
                record_entries{entry_index} = int_value;
            end

            record_index = digits_index(end)+1;
        case 'A' %string: Axxxxxxxx - always contains 8 characters
            STRING_LENGTH = 8;
            

            %%% entry per string 
            string_index = record_index + (1:STRING_LENGTH);

            %%% concatanate strings
            % test_index = record_index + 1 + STRING_LENGTH;
            % string_count = 1;
            % while data_record(test_index) == 'A' %concatenate adjacent strings
            %     string_count = string_count + 1;
            %     test_index = test_index + 1 + STRING_LENGTH;
            % end
            % string_index = zeros(1,STRING_LENGTH*string_count);
            % string_start_index = 1:9:(string_count*STRING_LENGTH);
            % for iString = 1:string_count
            %     string_span = (STRING_LENGTH*(iString-1)+1):(STRING_LENGTH*iString);
            %     string_index(string_span) = record_index + (string_start_index(iString):(string_start_index(iString) + STRING_LENGTH - 1));
            % end
            string_value = strip(data_record(string_index),"right",' ');

            entry_index = entry_index + 1;
            record_entries{entry_index} = string_value;

            record_index = string_index(end) + 1;
        case 'D' %double precision float: Dmx.xxxxxxxxxxxxxxxD±yy - where m is '-' if negative and y is the exponent 10^{±yy}
            FLOAT_LENGTH = 18;
            EXPONENT_LENGTH = 4;
            float_index = record_index + (1:FLOAT_LENGTH);
            float_coefficient = str2double(data_record(float_index));
            
            exponent_index = float_index(end) + (2:EXPONENT_LENGTH);
            float_exponent = str2double(data_record(exponent_index));

            float_value = float_coefficient*10^float_exponent;

            entry_index = entry_index + 1;
            record_entries{entry_index} = float_value;

            record_index = exponent_index(end) + 1;
        case ' '
            break
        otherwise
            %E for single precision float
            error("Unknown record type: '" + record_char + "'")

    end
end
% record_entries((entry_index + 1):end) = [];
end