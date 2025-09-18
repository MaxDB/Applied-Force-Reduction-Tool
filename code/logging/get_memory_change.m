function memory_change = get_memory_change
memory_data = get_free_memory;
memory_change = max(memory_data) - min(memory_data);
end