function [ gain ] = optical_LOS_gain_for_VLP( LED_location , receiver_location , parameter_m , detector_area)
% optical LOS gain for VLP

height_room = LED_location(3);
height_receiver = receiver_location(3);

distance = norm(receiver_location - LED_location);
gain = (parameter_m + 1)  * detector_area / (2 * pi * distance^2) * ( (height_room - height_receiver) / distance )^(parameter_m+1);

end

