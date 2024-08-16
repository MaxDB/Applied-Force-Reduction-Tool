function data = coco_energy_print(prob,data,command,LogLevel,varargin)
switch command
    case 'init'
        coco_print(prob,LogLevel,'   ENERGY')
    case 'data'
        bd_data = coco_bd_col(prob.bd, {"LAB","ENERGY","po.period"});
        label = bd_data(1,:);
        E = bd_data(2,:);
        
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
            num_points = size(label,2);
            if num_points == 1
                last_points = 1;
            else
                if label(end) > label(1)
                    last_points = num_points + [-1,0];
                else
                    last_points = [2,1];
                end
            end

            energy_percent = 100*E(last_points(end))/energy_limit;
            coco_print(prob, LogLevel, '% 9.1f',energy_percent)

            period = bd_data(3,:);
            coco_plot(period(last_points),E(last_points),energy_limit)
        end
        
end
end