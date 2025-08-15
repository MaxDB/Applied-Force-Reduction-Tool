function Model = model_eigenanalysis(Model,load_cache)
geometry_path = "geometry\" + Model.system_name + "\";
switch Model.Static_Options.static_solver
    case "abaqus"
        matrix_path = geometry_path + "matrices";
        matrices_loaded = isfile(matrix_path + ".mat") && load_cache;
        if matrices_loaded
            load(matrix_path,"M","K","matrix_bcs")

            logger("Matrices Loaded",3)
        else
            [M,K,matrix_bcs] = matrices_abaqus(Model.system_name);
            save(matrix_path,"M","K","matrix_bcs")
        end

    case "matlab"
        Analytic_Eom = load_analytic_system(geometry_path + Model.system_name);
        M = Analytic_Eom.linear_mass;
        K = Analytic_Eom.linear_stiffness;
        matrix_bcs = [];


end
Model.num_dof = size(M,1);

matrix_data = whos("M");
if matrix_data.bytes/1024 > Model.Static_Options.max_matrix_size
    Model.mass = Large_Matrix_Pointer(M,matrix_path,"M","save",0);
    Model.stiffness = Large_Matrix_Pointer(K,matrix_path,"K","save",0);
else
    Model.mass = M;
    Model.stiffness = K;
end


Model.dof_boundary_conditions = matrix_bcs;


if Model.Static_Options.load_custom_eigendata

    load("geometry\" + Model.system_name+"\custom_eigen_data.mat","custom_eval","custom_evec");

    r_modes = Model.reduced_modes;
    eval_r = custom_eval(r_modes);
    evec_r = custom_evec(:,r_modes);

else
    switch Model.Static_Options.additional_data
        case "perturbation"
            r_modes = Model.reduced_modes;
            L_modes = 1:Model.Static_Options.num_validation_modes;
            L_modes(ismember(L_modes,r_modes)) = [];
            Model.low_frequency_modes = L_modes;
            h_modes = [r_modes,L_modes];

            [evec,eval] = eigs(K,M,max(h_modes),"smallestabs");

            eVal_L = eval(L_modes,L_modes)*ones(length(L_modes),1);
            eVec_L = evec(:,L_modes);

            Model.low_frequency_eigenvalues = eVal_L;
            Model.low_frequency_eigenvectors = eVec_L;

        otherwise
            r_modes = Model.reduced_modes;

            [evec,eval] = eigs(K,M,max(r_modes),"smallestabs");


    end

    eval_r = eval(r_modes,r_modes)*ones(length(r_modes),1);
    evec_r = evec(:,r_modes);

end
Model.reduced_eigenvalues = eval_r;
matrix_data = whos("evec_r");
if matrix_data.bytes/1024 > Model.Static_Options.max_matrix_size
    data_path = Model.get_data_path;
    evec_r = Large_Matrix_Pointer(evec_r,data_path,"reduced_eigenvectors");
end
Model.reduced_eigenvectors = evec_r;


end