function coco_detected = coco_check()

try
coco_prob();
coco_detected = true;
catch exception
    coco_detected = false;
    if isequal(exception.identifier,'MATLAB:UndefinedFunction')
        fprintf(2,"Warning: Cannot detect COCO: ")
        fprintf('Please follow the <a href ="https://sourceforge.net/p/cocotools/wiki/installation/">installation guide</a>\n')
    else 
        rethrow(exception)
    end
end

if coco_detected
    fprintf('COCO detected\n')
end

end