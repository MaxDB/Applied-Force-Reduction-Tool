function [displacement,energy,additional_data,additional_data_time] = read_abaqus_fil_data(file_name,step_types,num_nodes,num_dimensions)

num_dofs = num_nodes*num_dimensions;
num_steps = length(step_types);
num_static_steps = sum(step_types(:,1) == "static");

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
if size(step_indicies,1) == 2*num_steps + 1
    step_indicies(1) = [];
else
    error("")
end
step_spans = zeros(num_steps,2);
step_spans(:,1) = step_indicies(1:2:(2*num_steps));
step_spans(:,2) = step_indicies(2:2:(2*num_steps));

for iStep = 1:num_steps
    step_span = step_spans(iStep,1):step_spans(iStep,2);
    step_type = step_types(iStep);
    switch step_type
        case "static"
            step_data = abaqus_data(step_span);
            energy_index = find_record_by_type(step_data,"energy",Model_Data);
            step_energy = read_energy(step_data(energy_index));
            
            
            disp_index = find_record_by_type(step_data,"displacement",Model_Data);
            step_displacement = read_displacements(step_data(disp_index));

            energy(iStep) = step_energy;
            displacement(:,iStep)  = reshape(step_displacement',num_dofs,1);
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
function displacement = read_displacements(disp_data)
    disp_table = split(disp_data,["D","I"+wildcardPattern(2)]);
    disp_mat = str2double(disp_table(:,4:end));
    %node_id = disp_mat(:,1); %for debugging
    
    disp_mat_length = size(disp_mat,2);   
    coefficient_index = 2:2:disp_mat_length;
    exponent_index = 3:2:disp_mat_length;

    displacement = disp_mat(:,coefficient_index).*10.^disp_mat(:,exponent_index);
end
%-------
function energy = read_energy(energy_data)
    energy_table = split(energy_data,["D","I"+wildcardPattern(2)]);
    energy_mat = str2double(energy_table(6:7));
    energy = energy_mat(1)*10^energy_mat(2);
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