function [x_grid,y_grid] = get_2d_plotting_data(Poly,Potential_Poly,energy_lim)


if nargin == 1 || isempty(Potential_Poly)
    GRID_DENSITY = 50;
    limits = Poly.input_limit;
    poly_bound = polyshape(limits');
    [x_lim,y_lim] = boundingbox(poly_bound);


    x = linspace(x_lim(1),x_lim(2),GRID_DENSITY);
    y = linspace(y_lim(1),y_lim(2),GRID_DENSITY);
    [X,Y] = meshgrid(x,y);

    X_BC = nan(GRID_DENSITY);
    Y_BC = nan(GRID_DENSITY);
    for iCol = 1:GRID_DENSITY
        x_vec = [X(:,iCol),Y(:,iCol)];
        valid_point = isinterior(poly_bound,x_vec);
        X_BC(valid_point,iCol) = X(valid_point,iCol);
        Y_BC(valid_point,iCol) = Y(valid_point,iCol);
    end

    x_grid = X_BC;
    y_grid = Y_BC;
else
    THETA_DENSITY = 100;
    RADIUS_DENSITY = 10;
    grid_density = ceil(sqrt(THETA_DENSITY*RADIUS_DENSITY));

    bounding_polygon = Poly.input_limit;
    r1_lim = [min(bounding_polygon(1,:)),max(bounding_polygon(1,:))];
    r2_lim = [min(bounding_polygon(2,:)),max(bounding_polygon(2,:))];

    r1_vec = linspace(r1_lim(1),r1_lim(2),grid_density);
    r2_vec = linspace(r2_lim(1),r2_lim(2),grid_density);
    [r1_grid,r2_grid] = meshgrid(r1_vec,r2_vec);

    r1_lin = reshape(r1_grid,1,grid_density^2);
    r2_lin = reshape(r2_grid,1,grid_density^2);
    V = Potential_Poly.evaluate_polynomial([r1_lin;r2_lin]);
    V = reshape(V,grid_density,grid_density);
    limit_contour = contourc(r1_vec,r2_vec,V,[energy_lim, energy_lim]);
    if floor(limit_contour(2,1)) ~= (size(limit_contour,2) - 1)
        warning("multiple energy limits")
    end
    limit_contour = limit_contour(:,2:(limit_contour(2,1)+1));
    limit_x = limit_contour(1,:);
    limit_y = limit_contour(2,:);
    % 

    limit_arc_length = estimate_arc_length(limit_x,limit_y);

    arc_length = linspace(0,limit_arc_length(end),THETA_DENSITY);
    arc_x = interp1(limit_arc_length,limit_x,arc_length);
    arc_y = interp1(limit_arc_length,limit_y,arc_length);

    arc_radius = sqrt(arc_x.^2 + arc_y.^2);
    arc_angle = atan2(arc_y,arc_x);

    radius_scale_factors = linspace(0,1,RADIUS_DENSITY);

    theta_grid = repmat(arc_angle',1,RADIUS_DENSITY);
    radius_grid = repmat(arc_radius',1,RADIUS_DENSITY).*radius_scale_factors;

    x_grid = radius_grid.*cos(theta_grid);
    y_grid = radius_grid.*sin(theta_grid);
end
end


function arc_length = estimate_arc_length(x,y)
diff_x = diff(x);
diff_y = diff(y);
incremental_length = [0,sqrt(diff_x.^2 + diff_y.^2)];
arc_length = cumsum(incremental_length);
end