function model_calibration_plot(mode,sep_id,iMode,r,f,f_limit,E,obj)
PLOT_LEVEL = 1;

load("data\plot_level.mat","plotting_level")
if plotting_level < PLOT_LEVEL
    return
end


figure
tiledlayout(1,2)
nexttile
box on
xlabel("r_" + mode)
ylabel("f_" + mode)

hold on
for iSep = 1:2
    sep_span = sep_id == (iSep+2*(iMode-1));
    plot([0,r(:,sep_span)],[0,f(:,sep_span)],'.-')
end

for iSep = 1:2
    plot(gca().XLim,f_limit(iSep)*[1,1],'k-')
end
plot(0,0,'k.','MarkerSize',10)
hold off

nexttile
box on
xlabel("r_" + mode)
ylabel("V")

hold on

for iSep = 1:2
    sep_span = sep_id == (iSep+2*(iMode-1));
    plot([0,r(:,sep_span)],[0,E(sep_span)],'.-')
end
plot(gca().XLim,obj.fitting_energy_limit*[1,1],'k-')
plot(0,0,'k.','MarkerSize',10)
hold off

end