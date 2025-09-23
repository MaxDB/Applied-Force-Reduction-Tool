function ax=plot_memory(varargin)

switch nargin
    case 0
        [memory_data,sample_delay] = get_free_memory;
    case 1
        memory_data = varargin{1};
        sample_delay = varargin{2};
end
time = 0:sample_delay:(length(memory_data)-1);


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