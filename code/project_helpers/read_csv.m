function output = read_csv(file_path)
% col delimeters: commas
% row delimiters: line break

ROW_DELIM = "\n";
COL_DELIM = ",";
size = [2,inf];

fID = fopen(file_path);
try 
    csv_data = fscanf(fID,"%f"+COL_DELIM+"%f",size);
    fclose(fID);
catch excemption
    fclose(fID);
    rethrow(excemption);
end

output=csv_data;

end