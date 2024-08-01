function set_logging_level(logging_level)
if ~exist("data","dir")
    mkdir("data")
end
save("data\log_level.mat","logging_level")
end