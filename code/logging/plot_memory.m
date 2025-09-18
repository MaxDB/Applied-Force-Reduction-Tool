function ax=plot_memory
SAMPLE_DEALY = 1;

memory_data = get_free_memory;
time = 0:SAMPLE_DEALY:(length(memory_data)-1);


figure
ax = axes;
plot(ax,time,memory_data,"-");
box(ax,"on")
xlabel(ax,"Time (s)")
ylabel(ax,"Free memory (GB)")
y_limit = ax.YLim;
y_limit(1) = 0;
ylim(ax,y_limit);


end