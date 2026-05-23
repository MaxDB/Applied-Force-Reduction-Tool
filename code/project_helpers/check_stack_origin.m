function is_live_script = check_stack_origin
    stack = dbstack;
    original_script = stack(end).file;
    is_live_script = startsWith(original_script,"LiveEditor");
end