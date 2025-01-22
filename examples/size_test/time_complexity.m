clear
% par_pool = gcp('nocreate');
% if isempty(par_pool)
%     parpool("Processes")
% end

%---------------------------------
% 0.01    --> 7,000
% 0.005   --> 24,500
% 0.003   --> 100,000
% 0.002   --> 300,000
% 0.00125 --> 1,000,000
% 0.001   --> 2,000,000
seed_sizes = linspace(0.01,0.00125,101);
% dof ≈ 0.0458 * seed_size ^ -2.53
%-------------------------------

num_seeds = length(seed_sizes);
num_dof = zeros(1,num_seeds);
static_time = zeros(2,num_seeds);
for iSeed = 1:num_seeds
    seed_size = seed_sizes(iSeed);

    %mesh arch with a particular seed size
    num_dof(iSeed) = create_mesh(seed_size);

    % %create static data
    % static_time_start = tic;
    % create_static_data(1,"none")
    % static_time(1,iSeed) = toc(static_time_start);
    % %create backbone
    % 
    % %validate backbone
    % 
    % %create static data
    % static_time_start = tic;
    % create_static_data([1,6],"none")
    % static_time(2,iSeed) = toc(static_time_start);
end
