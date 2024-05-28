function New_Opts = update_options(Default_Opts,Current_Opts,Input_Opts)

if isempty(Current_Opts)
    New_Opts = Default_Opts;
else
    New_Opts = Current_Opts;
end
opt_fields = fieldnames(Input_Opts);

for iField = 1:length(opt_fields)
    opt_field = opt_fields{iField,1};
    if ~isfield(New_Opts,opt_field)
        error("Unknown option: " + opt_field)
    end
    New_Opts.(opt_field) = Input_Opts.(opt_field);
end
end