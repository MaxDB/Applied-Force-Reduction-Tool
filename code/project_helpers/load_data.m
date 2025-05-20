function mat = load_data(obj)
if class(obj) == "Large_Matrix_Pointer"
    mat = obj.load;
    return
end

mat = obj;

end