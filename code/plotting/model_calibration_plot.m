function model_calibration_plot(mode,sep_id,iMode,r,f,f_limit,E,f_poly,v_poly,obj)
PLOT_LEVEL = 1;

load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end


figure
tiledlayout(1,2)
ax = nexttile;
box on
xlabel("r_" + mode)
ylabel("f_" + mode)



hold(ax,"on")
for iSep = 1:2
    tag = "sep_"+iSep;
    f_poly{iSep}.plot_polynomial("axes",ax,"tag",tag);
    
    sep_span = sep_id == (iSep+2*(iMode-1));
    plot(ax,r(:,sep_span),f(:,sep_span),'x',"Tag",tag)
end

for iSep = 1:2
    plot(ax,ax.XLim,f_limit(iSep)*[1,1],'k-')
end
plot(ax,0,0,'k.','MarkerSize',10)
hold(ax,"off")

ax = nexttile;
box on
xlabel("r_" + mode)
ylabel("V")



hold on
for iSep = 1:2
    tag = "sep_"+iSep;
    v_poly{iSep}.plot_polynomial("axes",ax,"tag",tag);
    
    sep_span = sep_id == (iSep+2*(iMode-1));
    plot(r(:,sep_span),E(sep_span),'x',"Tag",tag)
end
plot(gca().XLim,obj.fitting_energy_limit*[1,1],'k-')
plot(0,0,'k.','MarkerSize',10)
hold off

end