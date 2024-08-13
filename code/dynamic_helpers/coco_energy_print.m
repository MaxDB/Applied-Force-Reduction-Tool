function data = coco_energy_print(prob,data,command,LogLevel,varargin)
switch command
    case 'init'
        coco_print(prob,LogLevel,'   ENERGY')
    case 'data'
        bd_data = coco_bd_col(prob.bd, {"LAB","ENERGY","po.period"});
        
        E = bd_data(2,:);
        period = bd_data(3,:);
        for iField = 1:length(prob.efunc.events)
            if prob.efunc.events(iField).par == "ENERGY"
                energy_limit = prob.efunc.events(iField).vals;
                break
            end
        end
        if isempty(E)
            energy_percent = "---";
            coco_print(prob, LogLevel, '%9s',energy_percent)
        else
            energy_percent = 100*E(end)/energy_limit;
            coco_print(prob, LogLevel, '% 9.1f',energy_percent)
            coco_plot(period,E,energy_limit)
        end
        
end
end