function set_logging_level(logging_level)
if ~isfolder("data")
    mkdir("data")
end
    save("data\log_level.mat","logging_level")
end