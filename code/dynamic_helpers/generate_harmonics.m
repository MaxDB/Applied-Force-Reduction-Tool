function get_harmonics = generate_harmonics(harmonics_data,type)

if nargin == 1
    type = "base";
end
num_harmonics = length(harmonics_data);

get_harmonics = @(t,freq) 0;

switch type
    case "base"
        scale_factor = @(t,freq) 1;
    case "diff_time"
        scale_factor = @(t,freq) freq;
    case "diff_period"
        scale_factor = @(t,freq) -freq.^2.*t/(2*pi);
end

for iHarmonic = 1:num_harmonics
    harmonic_data = harmonics_data(iHarmonic);
    harmonic_num = iHarmonic - 1;
    if harmonic_data == 0
        continue
    end

    switch type
        case {"diff_time","diff_period"}
            harmonic_data = -harmonic_num*1j*harmonic_data;
    end

 

    get_harmonics = @(t,freq) get_harmonics(t,freq) + real(harmonic_data)*cos(harmonic_num*freq.*t) + imag(harmonic_data)*sin(harmonic_num*freq.*t);
end
get_harmonics = @(t,freq) scale_factor(t,freq).*get_harmonics(t,freq);

end


function test_plot(get_harmonics,freq)
period = 2*pi/freq;
t = linspace(0,period,100);
y = get_harmonics(t,freq);
figure
plot(t,y)
end