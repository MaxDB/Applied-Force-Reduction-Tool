function set_visualisation_level(plotting_level)
if ~exist("data","dir")
    mkdir("data")
end
save("data\plot_level.mat","plotting_level")
end