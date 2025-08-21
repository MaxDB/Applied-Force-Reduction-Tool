function adjacent_sep_ratios = get_adjacent_sep_ratios(used_sep_ratios)
% find the adjacent, sep_order SEPs to used_sep_ratios
num_sep_ratios = size(used_sep_ratios,2);
num_r_modes = size(used_sep_ratios,1);
if num_r_modes == 1
    adjacent_sep_ratios = [];
    return
end
num_angles = num_r_modes - 1;

phi = zeros(num_angles,num_sep_ratios);
adjacent_sep_ratios = zeros(num_r_modes,0);
for iSep = 1:num_sep_ratios
    sep_ratio = used_sep_ratios(:,iSep);
    phi(:,iSep) = cartesian_to_angles(sep_ratio);
    
    phi_adj = get_adjacent_angles(phi(:,iSep));
    num_adj_seps = size(phi_adj,2);
    for jSep = 1:num_adj_seps
        unit_sep_adj = angles_to_cartesian(phi_adj(:,jSep));
        adjacent_sep_ratios = [adjacent_sep_ratios,unit_sep_adj]; %#ok<AGROW>
    end
end
adjacent_sep_ratios = uniquetol(adjacent_sep_ratios', "ByRows",1)';

end

function phi = cartesian_to_angles(unit_sep_ratio)
num_r_modes = size(unit_sep_ratio,1);

phi = zeros(num_r_modes-1,1);
cartesian_squared = unit_sep_ratio.^2;
for iMode = 1:(num_r_modes-1)
    if all(isapprox(unit_sep_ratio(iMode:end),0))
        phi(iMode) = 0;
        continue
    end
    
    if iMode < (num_r_modes - 1)
        coord_sum = sum(cartesian_squared((iMode+1):end));
        tan_numerator = sqrt(coord_sum);
    elseif iMode == (num_r_modes - 1)
        tan_numerator = unit_sep_ratio(end);
    end
    phi(iMode) = atan2(tan_numerator,unit_sep_ratio(iMode));
end

end

%---
function unit_sep_ratio = angles_to_cartesian(phi)
num_angles = size(phi,1);
num_r_modes = num_angles + 1;

unit_sep_ratio = zeros(num_r_modes,1);
for iAngle = 1:num_angles
    sin_prod = 1;
    for jAngle = 1:(iAngle-1)
        sin_prod = sin_prod*sin(phi(jAngle));
    end
    unit_sep_ratio(iAngle) = cos(phi(iAngle))*sin_prod;
end

unit_sep_ratio(end) = prod(sin(phi));

end

%---

function phi_adj = get_adjacent_angles(phi)
num_angles = size(phi,1);

sep_order = get_sep_order(phi) + 1;
angle_inc = pi*2^(-sep_order);



increments = repelem({[-1,0,1]},1,num_angles);
incremental_combinations = table2array(combinations(increments{:}))';
num_combinations = size(incremental_combinations,2);
phi_adj = zeros(num_angles,num_combinations-1);
adj_counter = 0;
for iComb = 1:num_combinations
    inc_combo = incremental_combinations(:,iComb);
    if all(inc_combo == 0)
        continue
    end

    adj_counter = adj_counter + 1;
    phi_adj(:,adj_counter) = phi + inc_combo*angle_inc;
end

phi_adj = uniquetol(phi_adj', "ByRows",1)';

end

function order = get_sep_order(phi)
MAX_ORDER = 10;

for iOrder = 1:MAX_ORDER
    order = iOrder;
    base = pi*2^(-order);
    if isapprox(mod(phi/base,1),0,"tight") || isapprox(mod(phi/base,1),1,"tight")
        return
    end
end
if order == MAX_ORDER
    error("Could not detect SEP order")
end
end