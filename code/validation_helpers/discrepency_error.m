function error = discrepency_error(data_one,data_two)
    data_norm = max(abs(data_one),[],2);
    data_diff = abs(data_one-data_two);
    error = data_diff./data_norm;
end