function additional_disp = additional_physical_displacement(input_disp,Disp_Poly,Additional_Output,num_bc_dof,dof_bcs)
num_dof = num_bc_dof + numel(dof_bcs);
node_map = 1:num_dof;
node_map(dof_bcs) = [];
input_dimension = size(input_disp,1);
num_reduced_modes = Disp_Poly.input_dimension;

if isstring(Additional_Output.dof)
    if Additional_Output.dof == "all"
        fom_dof = 1:num_dof;        
        control_dof = 1:num_bc_dof;
    end
else
    fom_dof = Additional_Output.dof;
    control_dof = find(node_map == fom_dof);
end



switch input_dimension
    case num_reduced_modes  %input_disp = r
        phy_disp = Disp_Poly.evaluate_polynomial(input_disp,control_dof); 
    case num_dof            %input_disp = x
        phy_disp = input_disp(fom_dof,:);
    case num_bc_dof         %input_disp = x_bc
        phy_disp = input_disp(control_dof,:);
end

switch Additional_Output.type
    case "max"
        additional_disp = max(abs(phy_disp),[],2);
        if isstring(Additional_Output.dof)
            if Additional_Output.dof == "all"
                additional_disp = max(additional_disp);
            end
        end
    case "amplitude"
        additional_disp = (max(phy_disp,[],2) - min(phy_disp,[],2))/2;
end
end