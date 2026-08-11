function Analytic_Eom = load_analytic_system(system_path)
if ~isfile(system_path)
    if isfile(system_path + ".m")
        system_path = system_path + ".m";
    elseif isfile("geometry\" + system_path + "\" + system_path + ".m")
        system_path = "geometry\" + system_path + "\" + system_path + ".m";
    end
end
run(system_path)
if ~exist("Analytic_Eom","var")
    error("Analytic system isn't created. Ensure 'Analytic_Eom' is not renamed")
end
end