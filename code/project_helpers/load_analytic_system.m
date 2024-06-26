function Analytic_Eom = load_analytic_system(system_path)
    run(system_path)
    if ~exist("Analytic_Eom","var")
        error("Analytic system isn't created. Ensure 'Analytic_Eom' is not renamed")
    end
end