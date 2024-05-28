function data = coco_energy_print(prob,data,command,LogLevel,varargin)
switch command
    case 'init'
        coco_print(prob,LogLevel,'   ENERGY')
    case 'data'
        E = coco_bd_col(prob.bd, {"ENERGY"});
        for iField = 1:length(prob.efunc.events)
            if prob.efunc.events(iField).par == "ENERGY"
                energy_limit = prob.efunc.events(iField).vals;
                break
            end
        end
        coco_print(prob, LogLevel, '% 9.1f',100*E(end)/energy_limit)
end
end