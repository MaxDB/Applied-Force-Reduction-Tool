function ax = coco_plot(period,energy,energy_limit,ax)
PLOT_LEVEL = 1;
COLOUR = get_plot_colours(1);
MARKER_COLOUR = get_plot_colours(3);
LINE_WIDTH = 1.5;
MARKER_SIZE = 6;

plotting_level = load_log_level("plot");
if plotting_level < PLOT_LEVEL
    return
end

line_settings = {"Color",COLOUR,"LineWidth",LINE_WIDTH,"Marker",".","MarkerEdgeColor",MARKER_COLOUR};
marker_settings = {"Color",COLOUR,"LineWidth",LINE_WIDTH,"Marker","o","MarkerSize",MARKER_SIZE,"MarkerFaceColor",COLOUR,"MarkerEdgeColor",MARKER_COLOUR};

frequency = 2*pi./period;

% if isscalar(energy)
if isempty(ax)
    if check_stack_origin
        fig = figure('Internal',true);   % suppresses inline figure in live script
        set(fig,'Visible','on');   % off by default in live scripts
    else
        fig = figure;
    end
    fig.Name = "coco_plot";
    ax = axes(fig);
    hold(ax,"on")
    plot(ax,frequency,100*energy/energy_limit,marker_settings{:})
    % ylim([0,100])
    xlabel(ax,"Frequency (rad/s)")
    ylabel(ax,"Energy (% of max)")
    box(ax,"on")
    return
end

% figures = get(groot, 'Children');
% num_figures = size(figures,1);
% for iFig = 1:num_figures
%     fig = figures(iFig);
%     if fig.Name == "coco_plot"
%         ax = fig.Children;
%         break
%     end
% end
delete(ax.Children(1))

x_plot = frequency(1,(end-1):end);
y_plot = 100*energy(1,(end-1):end)/energy_limit;
plot(ax,x_plot,y_plot,line_settings{:})
plot(ax,x_plot(2),y_plot(2),marker_settings{:})
drawnow

end