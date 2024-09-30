function scaled_force_ratios = scale_sep_ratios(unit_force_ratios,force_magnitude)
if isempty(force_magnitude)
    force_magnitude = ones(size(unit_force_ratios));
end
num_modes = size(unit_force_ratios,1);

force_sign = sign(unit_force_ratios);

scaled_force_ratios = zeros(size(unit_force_ratios));
for iMode = 1:num_modes
    positive_index = force_sign(iMode,:) == 1;
    negative_index = force_sign(iMode,:) == -1;

    positive_force = unit_force_ratios(iMode,positive_index)*force_magnitude(iMode,1);
    negative_force = unit_force_ratios(iMode,negative_index)*abs(force_magnitude(iMode,2));

    scaled_force_ratios(iMode,positive_index) = positive_force;
    scaled_force_ratios(iMode,negative_index) = negative_force;
end

end