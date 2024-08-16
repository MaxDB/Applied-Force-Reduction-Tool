function validation_plot(plot_part,varargin)
%plot results of one iteration of validation alogirithm
PLOT_LEVEL = 3;
load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end



switch plot_part
    case 1 
        Static_Data = varargin{1,1};
        rom = varargin{1,2};

        num_modes = length(Static_Data.Model.reduced_modes);
        if num_modes > 2
            return
        end

        %plot coupling with worst fit
        plot_index = varargin{1,3};

        fig1 = figure;
        fig1_id = fig1.Number;
        tiledlayout("flow")
        ax = cell(1,2);
        ax{1,1} = nexttile;
        ax{1,2} = nexttile;
        
        % Static_Data.plot_condensed_displacement(plot_index,ax(1,1));
        ax = plot_static_data("displacement",Static_Data,"outputs",plot_index,"axes",ax{1,1});
        rom.Physical_Displacement_Polynomial.plot_polynomial(ax(1,1),plot_index);
        title(ax{1,1},"Degree = " + rom.Physical_Displacement_Polynomial.polynomial_degree);
        
        
        % Static_Data.plot_condensed_displacement(plot_index,ax{1,2});
        ax{1,2} = plot_static_data("displacement",Static_Data,"outputs",plot_index,"axes",ax{1,2});

        %plot energy and limit details
        
        plot_static_data("force-displacement",Static_Data);
        fig2 = gcf;
        fig2_id = fig2.Number;
        r_modes = Static_Data.Model.reduced_modes;
        num_modes = length(r_modes);
        energy_lim = Static_Data.Model.energy_limit;
        plotting_energy_lim = Static_Data.Model.fitting_energy_limit;

        hold on
        switch num_modes
            case 1
                
                xlabel("f_{" + r_modes(1) + "}")
            case 2
                xlabel("f_{" + r_modes(1) + "}")
                ylabel("f_{" + r_modes(2) + "}")

                x_lim = ax.XLim;
                y_lim = ax.YLim;
                [X,Y] = meshgrid(x_lim,y_lim);

                E_lim = zeros(size(X)) + energy_lim;
                plotting_E_lim = zeros(size(X)) + plotting_energy_lim;

                mesh(X,Y,E_lim,'FaceAlpha',0.1,'FaceColor','r')
                mesh(X,Y,plotting_E_lim,'FaceAlpha',0.1,'FaceColor','b')
        end
        hold off

        
        save("data\current_figure_info.mat","fig1_id","fig2_id","plot_index")
    case 2
        LINE_WIDTH = 1;

        r = varargin{1,1};
        num_modes = size(r,1);
        if num_modes > 2
            return
        end
        rom_one = varargin{1,2};
        rom_two = varargin{1,3};

        load("data\current_figure_info.mat","fig1_id","fig2_id","plot_index")
        figs = groot().Children;
        num_figs = length(figs);
        for iFig = 1:num_figs
            fig = figs(iFig);
            if fig.Number == fig1_id
                fig1 = fig;
            end
            if fig.Number == fig2_id
                fig2 = fig;
            end
        end

        
        

        T = fig1.Children;
        if length(T.Children) == 2
            %add legend
            ax1 = T.Children(2);
            ax2 = T.Children(1);
            hold(ax2,"on")
            L1 = plot(ax2,0,0,'k-',"LineWidth",LINE_WIDTH);
            L2 = plot(ax2,0,0,'r-.',"LineWidth",LINE_WIDTH);
            hold(ax2,"off")
            legend([L1,L2], ...
                ["Degree = " + rom_one.Condensed_Displacement_Polynomial.polynomial_degree,...
                "Degree = " + rom_two.Condensed_Displacement_Polynomial.polynomial_degree],...
                "Location","best","AutoUpdate","off")
        else
            ax1 = T.Children(3);
            ax2 = T.Children(2);
        end

        
        

        

        theta_one = rom_one.Condensed_Displacement_Polynomial.evaluate_polynomial(r,plot_index);
        theta_two = rom_two.Condensed_Displacement_Polynomial.evaluate_polynomial(r,plot_index);

        hold(ax1,"on")
        switch num_modes
            case 1
                plot(ax1,r(1,:),theta_one,'k-',"LineWidth",LINE_WIDTH)
            case 2
                plot3(ax1,r(1,:),r(2,:),theta_one,'k-',"LineWidth",LINE_WIDTH)
        end
        hold(ax1,"off")
        
        hold(ax2,"on")
        switch num_modes
            case 1
                plot(ax2,r(1,:),theta_one,'k-',"LineWidth",LINE_WIDTH)
                plot(ax2,r(1,:),theta_two,'r-.',"LineWidth",LINE_WIDTH)
            case 2
                plot3(ax2,r(1,:),r(2,:),theta_one,'k-',"LineWidth",LINE_WIDTH)
                plot3(ax2,r(1,:),r(2,:),theta_two,'r-.',"LineWidth",LINE_WIDTH)
        end
        hold(ax2,"off")
        
        %--------------------------------------------------------------------------%
        if isscalar(fig2.Children)
            %add legend
            ax3 = fig2.Children;
            hold(ax3,"on")
            L1 = plot(ax3,0,0,'k-',"LineWidth",LINE_WIDTH);
            L2 = plot(ax3,0,0,'r-.',"LineWidth",LINE_WIDTH);
            hold(ax3,"off")
            legend([L1,L2], ...
                ["Degree = " + rom_one.Force_Polynomial.polynomial_degree,...
                "Degree = " + rom_two.Force_Polynomial.polynomial_degree],...
                "Location","best","AutoUpdate","off")
        else
            ax3 = fig2.Children(2);
        end

        energy_one = rom_one.Potential_Polynomial.evaluate_polynomial(r);
        energy_two = rom_two.Potential_Polynomial.evaluate_polynomial(r);

        force_one = rom_one.Force_Polynomial.evaluate_polynomial(r);
        % force_two = rom_two.Force_Polynomial.evaluate_polynomial(r);

        hold(ax3,"on")
        switch num_modes
            case 1
                plot(ax3,force_one(1,:),energy_one,'k-',"LineWidth",LINE_WIDTH)
                plot(ax3,force_one(1,:),energy_two,'r-.',"LineWidth",LINE_WIDTH)
            case 2
                plot3(ax3,force_one(1,:),force_one(2,:),energy_one,'k-',"LineWidth",LINE_WIDTH)
                plot3(ax3,force_one(1,:),force_one(2,:),energy_two,'r-.',"LineWidth",LINE_WIDTH)
        end
        hold(ax3,"off")
end


end