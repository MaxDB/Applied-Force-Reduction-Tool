clear
data_path = "data\size_data";
%---------------------------------
% 0.01    --> 7,000
% 0.005   --> 24,500
% 0.0044 is too big
% 0.003   --> 100,000
% 0.002   --> 300,000
% 0.00125 --> 1,000,000
% 0.001   --> 2,000,000
% seed_sizes = [0.0044,0.003,0.00195,0.00125];
seed_sizes = [0.003,0.00195,0.00125];
% dof ≈ 0.0458 * seed_size ^ -2.53
%-------------------------------

%-- setup
if ~isfolder(data_path)
    mkdir(data_path)
end

par_pool = gcp('nocreate');
if isempty(par_pool)
    parpool("Processes")
end
%--
Size_Data.seed_sizes = seed_sizes;

num_seeds = length(seed_sizes);
num_dof = zeros(1,num_seeds);
static_time = zeros(2,num_seeds);
for iSeed = 1:num_seeds
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);

    %create static data
    static_time_start = tic;
    create_static_data(1,"none")
    static_time(1,iSeed) = toc(static_time_start);
    %create backbone

    %validate backbone

    %create static data
    static_time_start = tic;
    create_static_data([1,6],"none")
    static_time(2,iSeed) = toc(static_time_start);

    Size_Data.static_time = static_time;
    Size_Data.num_dof = num_dof;

    save(data_path + "\size_data","Size_Data") 

end
