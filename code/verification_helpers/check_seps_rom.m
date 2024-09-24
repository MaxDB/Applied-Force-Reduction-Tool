function [sep_force_error,sep_disp_error,sep_limit,error_modes] = check_seps_rom(rom_one,rom_two,force_ratios,validated_dofs)
num_seps = size(force_ratios,2);

sep_limit = zeros(size(force_ratios,1),num_seps);
sep_force_error = zeros(1,num_seps);
sep_disp_error = zeros(1,num_seps);
error_modes = zeros(1,num_seps);
for iSep = 1:num_seps
    [disp_sep_one,force_sep_one] = find_sep_rom(rom_one,force_ratios(:,iSep));
    sep_limit(:,iSep) = force_sep_one(:,end);
    
    energy_sep_one = rom_one.Potential_Polynomial.evaluate_polynomial(disp_sep_one);
    energy_limit = energy_sep_one > rom_one.Model.energy_limit;
    disp_sep_one(:,energy_limit) = [];
    force_sep_one(:,energy_limit) = [];

    condensed_sep_one = rom_one.Condensed_Displacement_Polynomial.evaluate_polynomial(disp_sep_one);
    condensed_sep_one = condensed_sep_one(validated_dofs,:);
    
    %check if agreement in rom_two
    force_sep_two = rom_two.Force_Polynomial.evaluate_polynomial(disp_sep_one);
    condensed_sep_two = rom_two.Condensed_Displacement_Polynomial.evaluate_polynomial(disp_sep_one);
    condensed_sep_two = condensed_sep_two(validated_dofs,:);

    force_error = discrepency_error(force_sep_one,force_sep_two);
    sep_force_error(1,iSep) = max(force_error,[],"all");
        
    condensed_error = max(discrepency_error(condensed_sep_one,condensed_sep_two),[],2);
    condensed_error(isnan(condensed_error)) = 0;
    [max_disp_error,max_disp_error_index] = max(condensed_error);
    sep_disp_error(1,iSep) = max_disp_error;
    error_modes(1,iSep) = validated_dofs(max_disp_error_index);

    validation_plot(2,disp_sep_one,rom_one,rom_two)
end
end