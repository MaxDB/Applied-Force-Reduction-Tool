function force_ratios = add_sep_ratios(r_num,index,found_force_ratios)

if r_num == 1
    force_ratios = [1,-1];
    if nargin == 3 && ~isempty(found_force_ratios)
        force_ratios = setdiff(force_ratios',found_force_ratios',"rows")';
    end
    return
end

phi_length = 2^index + 1;
theta_length = 2*(phi_length - 1);

phi = linspace(0,1,phi_length);
theta = linspace(0,2,theta_length+1);
theta(end) = [];

num_phi_coords = r_num - 2;
num_coords = num_phi_coords + 1;
num_combs = theta_length*phi_length^(num_phi_coords);
force_ratios = zeros(r_num,num_combs);

combination_input = cell(1,num_coords);

for iPhi_coord = 1:num_phi_coords
    combination_input{1,iPhi_coord} = phi;
end
combination_input{1,end} = theta;


coord_combs = table2array(combinations(combination_input{:}));


for iComb = 1:num_combs
    for iMode = 1:r_num
        r_temp = 1;
        for iSin_term = 1:(iMode-1)
            r_temp = r_temp*sinpi(coord_combs(iComb,iSin_term));
        end
        if iMode < r_num
            r_temp = r_temp*cospi(coord_combs(iComb,iMode));
        end
        force_ratios(iMode,iComb) = r_temp;
    end
end
force_ratios = round(force_ratios,5);

force_ratios = unique(force_ratios', 'rows')';

if nargin == 3 && ~isempty(found_force_ratios)
    force_ratios = setdiff(force_ratios',found_force_ratios',"rows")';
end
end

