function animate_whisker(Poly,Nc_Data,period,varargin)
num_args = length(varargin);
if mod(num_args,2) == 1
    error("Invalid keyword/argument pairs")
end
keyword_args = varargin(1:2:num_args);
keyword_values = varargin(2:2:num_args);

Potential_Poly = [];

for arg_counter = 1:num_args/2
    switch keyword_args{arg_counter}
        case "energy"
            [Potential_Poly,limit] = keyword_values{arg_counter}{:};
    end
end

%-----
fps = 24;
mark_limit = ~isempty(Potential_Poly);

num_modes = Poly.input_dimension;
num_applied_forces = Nc_Data.num_applied_forces;
num_r_modes = num_modes - num_applied_forces;
r_limit = Poly.input_limit(1:num_r_modes,:);

r_lim = [min(r_limit),max(r_limit)];
num_r_points = 100;
r = linspace(r_lim(1),r_lim(2),100);

y_limit = Poly.evaluate_polynomial(Poly.input_limit);
y_lim = [min(y_limit),max(y_limit)]*1.5;
if y_lim(1) > 0
    y_lim(1) = 0;
end


num_time_points = fps*period;
t = linspace(0,period,num_time_points+1);
t(end) = [];
tWave =   abs(mod(t-period/4,period)-period/2) - period/4;
pWave = 2*pi/period * tWave;


fig = figure;
ax = axes;
fig.Visible = "off";

frames(num_time_points) = struct('cdata',[],'colormap',[]);
Bar = Progress_Bar(num_time_points);
for iTime = 1:num_time_points
    p = pWave(iTime);
    poly_in = [r;ones(num_applied_forces,num_r_points)*p];
    poly_out = Poly.evaluate_polynomial(poly_in);
    if mark_limit
        potential = Potential_Poly.evaluate_polynomial(poly_in);
        over_limit = potential > limit;
    end
    if iTime == 1
        line = plot(ax,r,poly_out);
        if mark_limit
            hold(ax,"on")
            limit_markers = plot(ax,r(:,over_limit),poly_out(:,over_limit),"ro");
            hold(ax,"off")
        end
        time_text = annotation("textbox",[0.8,0.07,0.1,0.1],"String","t=0","EdgeColor","none");
        ylim(ax,y_lim)
        xlim(ax,ax.XLim)
        xlabel("r")
        title(ax,"Ctrl-C to exit")
     
    else
        line.YData = poly_out;
        limit_markers.XData = r(:,over_limit);
        limit_markers.YData = poly_out(:,over_limit);

        time_text.String = sprintf("t=%.1f/%.1f",t(iTime),period);
    end
    drawnow
    frames(iTime) = getframe(ax);
    Bar = Bar.increment_bar;
end
fig.Visible = "on";
movie(ax,frames,10,fps)
