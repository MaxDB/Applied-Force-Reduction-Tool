function sparse_matrix = load_mtx(file_path)
mtx_id = fopen(file_path);
mtx_data = textscan(mtx_id,'%f %f %f','delimiter','\n');
fclose(mtx_id);

sparse_matrix = [mtx_data{:}];
end