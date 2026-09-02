function [Validation_Rom, Verification_Rom] = prep_validation_roms(Static_Data)
validation_rom_start = tic;


degree_shift = 0;
Static_Data.Dynamic_Validation_Data.degree = Static_Data.Dynamic_Validation_Data.degree + degree_shift;
Static_Data.verified_degree(2) = Static_Data.verified_degree(2) + 2;
Validation_Rom = Reduced_System(Static_Data,"id",2);
Static_Data.Dynamic_Validation_Data.degree = Static_Data.Dynamic_Validation_Data.degree - degree_shift;


verification_rom_start = tic;
degree_shift = 2;
Static_Data.Dynamic_Validation_Data.degree = Static_Data.Dynamic_Validation_Data.degree + degree_shift;
Verification_Rom = Reduced_System(Static_Data,"id",3);
Static_Data.Dynamic_Validation_Data.degree = Static_Data.Dynamic_Validation_Data.degree - degree_shift;
verification_rom_time = toc(verification_rom_start);

log_message = sprintf("Verification ROM created: %.1f seconds" ,verification_rom_time);
logger(log_message,3)


validation_rom_time = toc(validation_rom_start);
log_message = sprintf("Validation ROM created: %.1f seconds" ,validation_rom_time);
logger(log_message,2)
end