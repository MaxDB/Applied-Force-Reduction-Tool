function Dyn_Data = one_mode_rom_validation(Dyn_Data,step)

switch step
    case 1
        compare_validation(Dyn_Data,"validation error",1,1:10);
    case 2
        compare_validation(Dyn_Data,"validation error",2,[3,6]);
end

end