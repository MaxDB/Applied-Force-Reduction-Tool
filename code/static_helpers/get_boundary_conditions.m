function boundary_conditions = get_boundary_conditions(geometry)
BOUNDARY_LINE = "*Boundary";
TERMINATION = "*";

num_geometry_lines = size(geometry,1);

boundary_lines = find(startsWith(geometry,BOUNDARY_LINE,"IgnoreCase",true));
num_lines = size(boundary_lines,1);
boundary_conditions = {};
for iLine = 1:num_lines
    boundary_line = boundary_lines(iLine);
    end_line = num_geometry_lines;
    for jLine = (boundary_line+1):end_line
        geometry_line = geometry{jLine};
        if startsWith(geometry_line,TERMINATION)
            end_line = jLine - 1;
            break
        end
    end
    boundary_span = (boundary_line+1):end_line;
    boundary_condition = [{convertStringsToChars(BOUNDARY_LINE + ", op=NEW")};geometry(boundary_span)];
    boundary_conditions = [boundary_conditions;boundary_condition]; %#ok<AGROW>
end

end