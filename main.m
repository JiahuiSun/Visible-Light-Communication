clear
clc
close all
%% parameters of signal
N_FFT = 128;
% bandwidth = 10*1e6; % full bandwith, unit: Hz;
bandwidth = 100*1e6; % full bandwith, unit: Hz;
K = 0:N_FFT-1;
N_CP = 16;
Num_frame = 2000;
M_QAM = 16;
run_times = 30;
Num_LED = 4;
PA_times = 25;
h_normalized = [1, 0.5, 0.2, zeros(1, N_FFT-3)];
range = [2,17,113,128;% range of four subcarrier blocks
        18,33,97, 112;
        34,49,81, 96;
        50,64,66, 80];
%% parameters of LED and receiver
detector_area = 1e-4; % unit: m^2
semiangle_at_half_power = pi/3;
parameter_m = -log(2)/log(cos(semiangle_at_half_power));
% Power_0 = 20; % power of optical signal, unit: W
Power_0 = 10;
I0 = 4;
slope_LED_PI = Power_0 / I0;
responsivity = 0.6; % responsivity of receiver
height_room = 5;
%% four LEDs positions %%%
x_LED_1 = 2; % unit: m
y_LED_1 = 2; % unit: m
z_LED_1 = height_room; % unit: m
position_LED_1 = [x_LED_1 , y_LED_1 , z_LED_1];

x_LED_2 = -2; % unit: m
y_LED_2 = 2; % unit: m
z_LED_2 = height_room; % unit: m
position_LED_2 = [x_LED_2 , y_LED_2 , z_LED_2];

x_LED_3 = -2; % unit: m
y_LED_3 = -2; % unit: m
z_LED_3 = height_room; % unit: m
position_LED_3 = [x_LED_3 , y_LED_3 , z_LED_3];

x_LED_4 = 2; % unit: m
y_LED_4 = -2; % unit: m
z_LED_4 = height_room; % unit: m
position_LED_4 = [x_LED_4 , y_LED_4 , z_LED_4];

position_LED = [position_LED_1 ; position_LED_2 ; position_LED_3 ; position_LED_4];
position_receiver = [-4, -3, 1];
%% preallocation
distance = zeros(Num_LED, 1);
channel_optical_gain = zeros(Num_LED, 1);
h_actual = cell(Num_LED, 1);
rx_location_est_method_0 = zeros(2, run_times);
Delta_p_square = zeros(run_times, 1);
Delta_Di_square = zeros(3, run_times);
rx_location_est_method_0_PA = zeros(2, run_times);
Delta_p_square_PA = zeros(run_times, 1);
Delta_Di_square_PA = zeros(3, run_times);
rx_location_est_method_0_3LED = zeros(2, run_times);
rx_location_est_method_0_PA_tmp = zeros(2, PA_times);
Delta_p_square_3LED = zeros(run_times, 1);
Delta_Di_square_3LED = zeros(3, run_times);
%% set up channel
for i_LED = 1:Num_LED
    distance(i_LED) = norm(position_receiver - position_LED(i_LED,:));
    channel_optical_gain(i_LED) = optical_LOS_gain_for_VLP( position_LED(i_LED,:) , position_receiver , parameter_m , detector_area);
    h_actual{i_LED} = responsivity * slope_LED_PI * channel_optical_gain(i_LED) * h_normalized; %actual channel inpulse response
end
%% generate signal
for run = 1:run_times
    %% overdetermined VLP 
    [sig_ofdm_time, sig_ofdm_freq] = ofdmMod(Num_frame, N_CP, M_QAM, Num_LED, N_FFT, range, h_actual);
    noise_power = noise_p(channel_optical_gain, Power_0, bandwidth);
    [noise_row , noise_column] = size(sig_ofdm_time{1});
    noise = 1/sqrt(2) * sqrt(noise_power) * (randn(noise_row , noise_column) + 1j * randn(noise_row , noise_column));
    signal_received_with_noise =  (sig_ofdm_time{1} + sig_ofdm_time{2} + sig_ofdm_time{3} + sig_ofdm_time{4})+ noise;
    % signal_received_with_noise =  (sig_ofdm_time{1} + sig_ofdm_time{2} + sig_ofdm_time{3} + sig_ofdm_time{4});

    recovered_signal = 1 / sqrt(N_FFT) * fft(signal_received_with_noise(:, N_CP+1:end), N_FFT, 2);
    % scatterplot(recovered_signal(:, 18));
    tx_power = mean(abs(sig_ofdm_freq{1}+sig_ofdm_freq{2}+sig_ofdm_freq{3}+sig_ofdm_freq{4}).^2, 1);
    rx_power = mean(abs(recovered_signal).^2, 1);

%     SNR_dB = 10*log10(rx_power./noise_power);
%     figure();
%     stem(K, SNR_dB); hold on; grid on;
%     xlabel('K'); ylabel('SNR/dB');
%     title('SNR before PA');

    H_frequency_gain_square = rx_power ./ tx_power / (responsivity * slope_LED_PI)^2;

%     figure();
%     stem(K, H_frequency_gain_square); grid on;
%     title('|HK|^2 before PA');

    [~, estimated_distance] = curveFit(Num_LED, K, range, H_frequency_gain_square);
    estimated_distance = abs(estimated_distance);
    [M, I0] = sort(estimated_distance);
    
    rx_location_est_method_0(:,run) = VLP_LSE_3LEDs( position_LED( I0(1) , :) , position_LED( I0(2) , :) , position_LED( I0(3) , :) ,...
                                                       estimated_distance(I0(1)) , estimated_distance(I0(2)) , estimated_distance(I0(3)));

    Delta_p_square(run) = norm([-4 ; -3] - rx_location_est_method_0(:,run)).^2;
    
    Delta_Di_square(:, run) = [estimated_distance(I0(1));estimated_distance(I0(2));estimated_distance(I0(3))].^2 - [distance(I0(1));distance(I0(2));distance(I0(3))].^2;
    %% determined VLP
    range_3LED = [2 22 108 128;% range of three subcarrier blocks
                  23 43 87 107;
                  44 64 66 86];
    position_3LED = [position_LED(I0(3),:); position_LED(I0(2),:); position_LED(I0(1),:)];
    h_actual_3LED = {h_actual{I0(3)}; h_actual{I0(2)}; h_actual{I0(1)}};
    Num_LED_3LED = 3;
    channel_optical_gain_3LED = [channel_optical_gain(I0(3)), channel_optical_gain(I0(2)), channel_optical_gain(I0(1))];
    
    [sig_ofdm_time, sig_ofdm_freq] = ofdmMod(Num_frame, N_CP, M_QAM, Num_LED_3LED, N_FFT, range_3LED, h_actual_3LED);
    noise_power_3LED = noise_p(channel_optical_gain_3LED, Power_0, bandwidth);
    [noise_row , noise_column] = size(sig_ofdm_time{1});
    noise_3LED = 1/sqrt(2) * sqrt(noise_power_3LED) * (randn(noise_row , noise_column) + 1j * randn(noise_row , noise_column));
    signal_received_with_noise_3LED =  (sig_ofdm_time{1} + sig_ofdm_time{2} + sig_ofdm_time{3})+ noise_3LED;
    % signal_received_with_noise =  (sig_ofdm_time{1} + sig_ofdm_time{2} + sig_ofdm_time{3});

    recovered_signal = 1 / sqrt(N_FFT) * fft(signal_received_with_noise_3LED(:, N_CP+1:end), N_FFT, 2);
    % scatterplot(recovered_signal(:, 18));
    tx_power = mean(abs(sig_ofdm_freq{1}+sig_ofdm_freq{2}+sig_ofdm_freq{3}).^2, 1);
    rx_power = mean(abs(recovered_signal).^2, 1);

%     SNR_dB = 10*log10(rx_power./noise_power_3LED);
%     figure();
%     stem(K, SNR_dB); hold on; grid on;
%     xlabel('K'); ylabel('SNR/dB');
%     title('SNR before PA');

    H_frequency_gain_square_3LED = rx_power ./ tx_power / (responsivity * slope_LED_PI)^2;
    
%     figure();
%     stem(K, H_frequency_gain_square_3LED); grid on;
%     title('|HK|^2 before PA');

    [estimated_gain, estimated_distance] = curveFit(Num_LED_3LED, K, range_3LED, H_frequency_gain_square_3LED);
    estimated_distance = abs(estimated_distance);
    [M_3LED, I0_3LED] = sort(estimated_distance);
    
    rx_location_est_method_0_3LED(:,run) = VLP_LSE_3LEDs( position_LED( I0(1), :) , position_LED( I0(2) , :) , position_LED( I0(3) , :) ,...
                                               estimated_distance(I0_3LED(1)) , estimated_distance(I0_3LED(2)) , estimated_distance(I0_3LED(3)));

    Delta_p_square_3LED(run) = norm([-4 ; -3] - rx_location_est_method_0_3LED(:,run)).^2;
    
    Delta_Di_square_3LED(:, run) = [estimated_distance(I0_3LED(1)); estimated_distance(I0_3LED(2)); estimated_distance(I0_3LED(3))].^2 ...
                                               - [distance(I0(1)); distance(I0(2)); distance(I0(3))].^2;    
    %% PA
    for n = 1:PA_times
        if n == 1
            [sig_ofdm_time_PA, sig_ofdm_freq_PA, a_square] = PA_method(Num_frame, N_CP, M_QAM, Num_LED_3LED, N_FFT, range_3LED, h_actual_3LED, H_frequency_gain_square_3LED);
        else
            [sig_ofdm_time_PA, sig_ofdm_freq_PA, a_square] = PA_method(Num_frame, N_CP, M_QAM, Num_LED_3LED, N_FFT, range_3LED, h_actual_3LED, H_frequency_gain_square_PA);            
        end
%         a_square_all{n} = a_square;
%         if n > 1
%             detla_a_square = a_square_all{n}-a_square_all{n-1};
%             figure();
%             plot(K, detla_a_square, '-r*'); hold on; grid on;
%             title('a square');
%         end
        
        [noise_row , noise_column] = size(sig_ofdm_time_PA{1});
        noise = 1/sqrt(2) * sqrt(noise_power_3LED) * (randn(noise_row , noise_column) + 1j * randn(noise_row , noise_column));
        signal_received_with_noise_PA =  (sig_ofdm_time_PA{1} + sig_ofdm_time_PA{2} + sig_ofdm_time_PA{3})+ noise;
        % signal_received_with_noise =  (sig_ofdm_time{1} + sig_ofdm_time{2} + sig_ofdm_time{3});

        recovered_signal_PA = 1 / sqrt(N_FFT) * fft(signal_received_with_noise_PA(:, N_CP+1:end), N_FFT, 2);
        % scatterplot(recovered_signal(:, 18));
        tx_power = mean(abs(sig_ofdm_freq_PA{1}+sig_ofdm_freq_PA{2}+sig_ofdm_freq_PA{3}).^2, 1);
        rx_power = mean(abs(recovered_signal_PA).^2, 1);

%         SNR_dB = 10*log10(rx_power./noise_power_3LED);
%         figure();
%         stem(K, SNR_dB); hold on; grid on;
%         xlabel('K'); ylabel('SNR/dB');
%         title('SNR after PA');

        H_frequency_gain_square_PA = rx_power ./ tx_power / (responsivity * slope_LED_PI)^2 ./ a_square;

%         figure();
%         stem(K, H_frequency_gain_square_PA); grid on;
%         title('|HK|^2 after PA');

        [estimated_gain_PA, estimated_distance_PA] = curveFit(Num_LED_3LED, K, range_3LED, H_frequency_gain_square_PA);
        estimated_distance_PA = abs(estimated_distance_PA);
        [M_PA, I0_PA] = sort(estimated_distance_PA);
        rx_location_est_method_0_PA_tmp(:, n) = VLP_LSE_3LEDs( position_LED(I0(1), :), position_LED(I0(2), :), position_LED(I0(3), :) ,...
                                                       estimated_distance_PA(I0_PA(1)), estimated_distance_PA(I0_PA(2)), estimated_distance_PA(I0_PA(3)));        
    end
%     rx_location_est_method_0_PA(:, run) = VLP_LSE_3LEDs( position_LED(I0(1), :), position_LED(I0(2), :), position_LED(I0(3), :) ,...
%                                                        estimated_distance_PA(I0_PA(1)), estimated_distance_PA(I0_PA(2)), estimated_distance_PA(I0_PA(3)));
    rx_location_est_method_0_PA(:, run) = mean(rx_location_est_method_0_PA_tmp(:, 3:PA_times), 2);                                         
    Delta_p_square_PA(run) = norm([-4 ; -3] - rx_location_est_method_0_PA(:,run)).^2;
    
    Delta_Di_square_PA(:, run) = [estimated_distance_PA(I0_PA(1)); estimated_distance_PA(I0_PA(2)); estimated_distance_PA(I0_PA(3))].^2 - [distance(I0(1));distance(I0(2));distance(I0(3))].^2;
end
%% over determined VLP error analyse
rx_pos_mean = mean(rx_location_est_method_0, 2);
E_Delta_p = norm( [-4;-3] - rx_pos_mean );
MSE_pos = mean(Delta_p_square);

figure();
plot(rx_location_est_method_0(1,:), rx_location_est_method_0(2,:), '.');
xlabel('x'); ylabel('y');
axis([-6 6 -6 6]);
hold on; grid on;
plot(rx_pos_mean(1), rx_pos_mean(2), 'r*');
legend('4led estimated positon', '4led mean position');

figure()
plot(1:run_times, Delta_p_square, '*'); 
xlabel('n'); ylabel('MSE:detla p');
hold on; grid on;
plot([0 run_times], [MSE_pos MSE_pos], 'r'); hold on;
legend('4led estimated error', '4led mean error');

%% determined and PA data analyse
rx_pos_mean_3LED = mean(rx_location_est_method_0_3LED, 2);
E_Delta_p_3LED = norm( [-4;-3] - rx_pos_mean_3LED );
MSE_pos_3LED = mean(Delta_p_square_3LED);

figure();
plot(rx_location_est_method_0_3LED(1,:), rx_location_est_method_0_3LED(2,:), '.');
xlabel('x'); ylabel('y');
axis([-6 6 -6 6]);
hold on; grid on;
plot(rx_pos_mean_3LED(1), rx_pos_mean_3LED(2), 'r*');
legend('3led estimated positon', '3led mean position');

figure()
plot(1:run_times, Delta_p_square_3LED, '*'); 
xlabel('n'); ylabel('error');
hold on; grid on;
plot([0 run_times], [MSE_pos_3LED MSE_pos_3LED], 'r'); hold on;
legend('3led estimated error', '3led mean error');

%% determined VLP data analyse
rx_position_mean_PA = mean(rx_location_est_method_0_PA, 2);
E_Delta_p_PA = norm( [-4;-3] - rx_position_mean_PA );
MSE_pos_PA = mean(Delta_p_square_PA);

figure();
plot(rx_location_est_method_0_PA(1,:), rx_location_est_method_0_PA(2,:), '.');
xlabel('x'); ylabel('y');
axis([-6 6 -6 6]);
hold on; grid on;
plot(rx_position_mean_PA(1), rx_position_mean_PA(2), 'r*');
legend('3led PA estimated positon', '3led PA mean position');

figure()
plot(1:run_times, Delta_p_square_PA, '*'); 
xlabel('n'); ylabel('error');
hold on; grid on;
plot([0 run_times], [MSE_pos_PA MSE_pos_PA], 'r'); hold on;
legend('3led PA estimated error', '3led PA mean error');

%% Delta_Di_square analyse
Delta_Di_square_mean = mean(Delta_Di_square, 2);
Delta_Di_square_mean_PA = mean(Delta_Di_square_PA, 2);
Delta_Di_square_mean_3LED = mean(Delta_Di_square_3LED, 2);
figure()
plot([1, 2, 3], Delta_Di_square_mean, '-ks');
xlabel('n'); ylabel('Delta Di square mean');
hold on; grid on;
plot([1 2 3], Delta_Di_square_mean_PA, '-ro');
hold on; grid on;
plot([1 2 3], Delta_Di_square_mean_3LED, '-gd');
legend('Delta Di square mean', 'Delta Di square mean PA', 'Delta Di square mean 3LED');