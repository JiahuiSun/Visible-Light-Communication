function [ location ] = VLP_LSE_3LEDs( LED_1_location , LED_2_location , LED_3_location , distance_1_estimated , distance_2_estimated , distance_3_estimated)
% location estimation due to least square estimation with 3 LEDs

A = [  2*(LED_2_location(1) - LED_1_location(1)) , 2*(LED_2_location(2) - LED_1_location(2)) 
         2*(LED_3_location(1) - LED_1_location(1)) , 2*(LED_3_location(2) - LED_1_location(2)) ];

B = [  (LED_2_location(1)^2 - LED_1_location(1)^2) + (LED_2_location(2)^2 - LED_1_location(2)^2) + (distance_1_estimated^2 - distance_2_estimated^2)
          (LED_3_location(1)^2 - LED_1_location(1)^2) + (LED_3_location(2)^2 - LED_1_location(2)^2) + (distance_1_estimated^2 - distance_3_estimated^2) ];
      
location = (A' * A) ^(-1) * A' * B;
% location(1) for x
% location(2) for y

end

