function verification_plot(Static_Data,Model,stage)
%plot results of one iteration of verification alogirithm
PLOT_LEVEL = 3;
FIG_NAME = "verification plot";
load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end

figs = groot().Children;
num_figs = length(figs);
verification_fig = [];
for iFig = 1:num_figs
    fig = figs(iFig);
    if fig.Name == FIG_NAME
        verification_fig = fig;
    end
end

if isempty(verification_fig)
    verification_fig = figure;
    verification_fig.Name = FIG_NAME;
    box on
end
set(groot(),"CurrentFigure",verification_fig)
ax = gca;
%---------------------

num_modes = size(Model.reduced_modes,2);
energy_lim = Model.energy_limit;
colour = get_plot_colours(stage);

switch num_modes
    case 1
        switch stage
            case 0
                plot_static_data("energy",Static_Data,"axes",ax,"plot seps",0,"colour",colour);
                hold(ax,"on")

                x_lim = ax.XLim;
                E_lim = zeros(size(x_lim)) + energy_lim;
                plot(x_lim,E_lim,'r-');
               

            otherwise
                hold(ax,"on")
                r = Static_Data.r;
                E = Static_Data.E;
                plot(r,E,".","Color",colour)
        end
    case 2
        switch stage
            case 0
                plot_static_data("energy",Static_Data,"axes",ax,"plot seps",0,"colour",colour);
                hold(ax,"on")

                x_lim = ax.XLim;
                y_lim = ax.YLim;
                [X,Y] = meshgrid(x_lim,y_lim);
                E_lim = zeros(size(X)) + energy_lim;
                mesh(X,Y,E_lim,'FaceAlpha',0.1,'FaceColor',"none");
               

            otherwise
                hold(ax,"on")
                r = Static_Data.r;
                E = Static_Data.E;
                plot3(r(1,:),r(2,:),E,".","Color",colour)
        end
end
drawnow
hold(ax,"off")
end