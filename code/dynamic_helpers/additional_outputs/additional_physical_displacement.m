function additional_disp = additional_physical_displacement(input_disp,Disp_Poly,Additional_Output,node_map)
control_dof = node_map(node_map(:,1) == Additional_Output.dof,2);

num_dof = node_map(end,1);
num_bc_dof = node_map(end,2);

input_dimension = size(input_disp,1);
num_reduced_modes = Disp_Poly.input_dimension;
switch input_dimension
    case num_reduced_modes  %input_disp = r
        phy_disp = Disp_Poly.evaluate_polynomial(input_disp,control_dof); 
    case num_dof            %input_disp = x
        phy_disp = input_disp(Additional_Output.dof,:);
    case num_bc_dof         %input_disp = x_bc
        phy_disp = input_disp(control_dof,:);
end

switch Additional_Output.type
    case "max"
        additional_disp = max(abs(phy_disp));
    case "amplitude"
        additional_disp = (max(phy_disp) - min(phy_disp))/2;
end
end