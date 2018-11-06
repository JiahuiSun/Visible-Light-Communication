function [noise_power] = noise_p(channel_optical_gain, Power_0, bandwidth)
%% parameters
    detector_area = 1e-4; % unit: m^2
    responsivity = 0.6; % responsivity of receiver
    electronic_charge = 1.6021892*1e-19; % unit: C
    current_2 = 0.562;
    Boltzman_constant = 1.3806505*1e-23; % unit: J/K
    capacitance_per_area = 112*1e-12/1e-4; %unit:F/m^2
    open_loop_voltage_gain = 10;
    Temperature = 298; % unit: K
    noise_factor = 1.5;
    FET_transconductance = 30*1e-3; % unit: S
    current_3 = 0.0868;
    
    received_optical_power =  channel_optical_gain * Power_0;
    total_received_optical_power = sum(received_optical_power); 
    thermal_noise =  8*pi*Boltzman_constant*Temperature*capacitance_per_area*detector_area*current_2*bandwidth^2/open_loop_voltage_gain + ...
        16*pi^2*Boltzman_constant*Temperature*noise_factor*capacitance_per_area^2*detector_area^2*current_3*(bandwidth)^3/FET_transconductance;
    shot_noise = 2 * electronic_charge * responsivity *  total_received_optical_power * bandwidth;
    noise_power = shot_noise + thermal_noise;
end