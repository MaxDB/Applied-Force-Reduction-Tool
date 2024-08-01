function set_visualisation_level(plotting_level)
if ~isfolder("data")
    mkdir("data")
end
save("data\plot_level.mat","plotting_level")
end