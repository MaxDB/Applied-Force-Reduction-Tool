function colours = get_plot_colours(colour_range)

% Set default plotting colours
% https://colorbrewer2.org/ -> qualitative Set3
COLOURS = [
0.552941176470588, 0.827450980392157, 0.780392156862745;
0.745098039215686, 0.729411764705882, 0.854901960784314;
0.984313725490196, 0.501960784313726, 0.447058823529412;
0.501960784313726, 0.694117647058824, 0.827450980392157;
0.992156862745098, 0.705882352941177, 0.384313725490196;
0.701960784313725, 0.870588235294118, 0.411764705882353;
0.988235294117647, 0.803921568627451, 0.898039215686275;
0.737254901960784, 0.501960784313726, 0.741176470588235;
0.8,               0.921568627450980, 0.772549019607843;
1,                 0.929411764705882, 0.435294117647059
];

COLOUR_MAP = [
    255,247,251;
    236,231,242;
    208,209,230;
    166,189,219;
    116,169,207;
    54,144,192;
    5,112,176;
    4,90,141;
    2,56,88
    ]/256;

BLACK = [0,0,0];
GREY = [217,217,217]/256;

num_colours = size(COLOURS,1);

if nargin == 0
    colours = COLOURS;
    return
end

if isstring(colour_range)
    if any(colour_range == ["grey","gray"])
        colours = GREY;
        return
    end

    if colour_range == "map"
        colours = COLOUR_MAP;
        return
    end
end

if colour_range == 0
    colours = BLACK;
    return
end

colours = COLOURS(mod(colour_range-1,num_colours)+1,:);
end