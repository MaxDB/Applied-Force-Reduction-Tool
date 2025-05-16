function reset_temp_directory
try
    rmdir('temp','s')
    mkdir('temp')
catch
    try
        pause(0.5)
        rmdir('temp','s')
        mkdir('temp')
    catch
        warning("Cannot delete all files in temp directory")
    end
end

% if ~isfolder("sparse_temp")
%     return
% end
% 
% try
%     rmdir('sparse_temp','s')
% catch
%     try
%         pause(0.5)
%         rmdir('sparse_temp','s')
%     catch
%         warning("Cannot delete all files in sparse temp directory")
%     end
% end

% if isfolder('temp')
%     iLoop = 1;
%     while iLoop < 100
%         try
%             rmdir('temp','s')
%             break
%         catch
%             if iLoop == 1
%                 fclose("all");
%             end
%             iLoop = iLoop + 1;
%             if iLoop > 10
%                 error("Can't delete temp directory")
%             end
%             pause(0.1)
%         end
%     end
% end
% mkdir('temp')
end