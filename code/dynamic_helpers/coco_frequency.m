function [data,frequency] = coco_frequency(prob,data,u)
period = u(end);
frequency = 2*pi/period;
end