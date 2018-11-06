function [sig_ofdm_time, sig_ofdm_freq, a_square] = PA_method(Num_frame, N_CP, M_QAM, Num_LED, N_FFT, range, h_actual, H_frequency_gain_square)
    K_QAM = log2(M_QAM);
    N_FFT_half = N_FFT/2;
    numBit = K_QAM*Num_frame*(N_FFT_half-1)/2;

    sig_ofdm_freq = cell(Num_LED, 1);
    signal_LED_time_no_CP = cell(Num_LED, 1);
    signal_LED_time_add_CP = cell(Num_LED, 1);
    signal_LED_time_serial_output = cell(Num_LED, 1);
    signal_LED_time_reflection = cell(Num_LED, 1);
    signal_LED_time_parallel_output = cell(Num_LED, 1);
    sig_ofdm_freq_mul_a = cell(Num_LED, 1);
    sig_ofdm_time = cell(Num_LED, 1);

    if M_QAM == 8
        Energy_M_QAM = 6;
    else
        if sqrt(M_QAM) == fix(sqrt(M_QAM))
            Energy_M_QAM = 2 / 3 * (M_QAM - 1); % Normalized coefficient
        else
            P_QAM = (log2(M_QAM /4)- 1)/2;
            P_QAM0 = (2^(P_QAM-1));
            P_QAM3 = 3 * P_QAM0;
            P_QAM2 = P_QAM3 - P_QAM0 ; % Intermediate variables
            Energy_M_QAM = ((2 * P_QAM3) * (P_QAM3 * (4 * P_QAM3^2-1)/3) ...
                - (2 * P_QAM0) * ( P_QAM3 * (4 * P_QAM3^2-1)/3 - P_QAM2 * (4 * P_QAM2^2-1)/3)) / (P_QAM3^2 - P_QAM0^2);
        end
    end

    BitSource = randi([0 , 1] , 1 , numBit);  % source bit
    BitSource_PPM = zeros(1 , 2 * numBit);     % source bit after PPM
    for jj=1:numBit;
        if(BitSource(jj)==0)
            BitSource_PPM(2*jj)=1;
        else
            BitSource_PPM(2*jj-1)=1;
        end
    end
    % QAM modulation with normalized power
    symbol_decimal = bi2de(reshape(BitSource_PPM , numel(BitSource_PPM)/K_QAM , K_QAM), 'left-msb');
    symbol_input = reshape(symbol_decimal, [], Num_frame)';
    symbol_MQAM_matrix = 1 / sqrt(Energy_M_QAM) * qammod(symbol_input, M_QAM, 0);
    a_square = zeros(1,N_FFT);
    a_power_all = (range(end,2)-range(1,1)+1)*2;

    for ii = range(1,1):range(end,2)
        a_square(ii) = a_power_all*(1/H_frequency_gain_square(ii)) / ...
            sum(1./[H_frequency_gain_square(range(1,1):range(end,2)), H_frequency_gain_square(range(end,3):range(1,4))]);
    end
    for ii = range(end,3):range(1,4)
        a_square(ii) = a_power_all*(1/H_frequency_gain_square(ii)) / ...
            sum(1./[H_frequency_gain_square(range(1,1):range(end,2)), H_frequency_gain_square(range(end,3):range(1,4))]);
    end
    
    for i_LED = 1:Num_LED 
        sig_ofdm_freq{i_LED} = zeros(Num_frame, N_FFT);
        sig_ofdm_freq_mul_a{i_LED} = zeros(Num_frame, N_FFT);
        sig_ofdm_freq{i_LED}(: , range(i_LED, 1) : range(i_LED, 2)) = symbol_MQAM_matrix(:, range(i_LED, 1)-1:range(i_LED, 2)-1);
        sig_ofdm_freq{i_LED}(: , range(i_LED, 3) : range(i_LED, 4)) = conj(sig_ofdm_freq{i_LED}(:, range(i_LED, 2):-1:range(i_LED, 1))); 
        
        for dd = range(i_LED, 1) : range(i_LED, 2)
            sig_ofdm_freq_mul_a{i_LED}(:,dd) = sqrt(a_square(dd)) * sig_ofdm_freq{i_LED}(:,dd);
        end
        for dd = range(i_LED, 3) : range(i_LED, 4)
            sig_ofdm_freq_mul_a{i_LED}(:,dd) = sqrt(a_square(dd)) * sig_ofdm_freq{i_LED}(:,dd);
        end
        
        signal_LED_time_no_CP{i_LED} = sqrt(N_FFT) * ifft(sig_ofdm_freq_mul_a{i_LED}, N_FFT, 2);
        signal_LED_time_add_CP{i_LED} = [signal_LED_time_no_CP{i_LED}(:, end+1-N_CP:end), signal_LED_time_no_CP{i_LED}];

        signal_LED_time_serial_output{i_LED} = reshape(signal_LED_time_add_CP{i_LED}', 1, []);                
        % signal transmission with reflection
        signal_LED_time_reflection{i_LED} = h_actual{i_LED}(1) * signal_LED_time_serial_output{i_LED} + ...
                                            h_actual{i_LED}(2) * [0 , signal_LED_time_serial_output{i_LED}(1:end-1)]+ ...
                                            h_actual{i_LED}(3) * [0 , 0 , signal_LED_time_serial_output{i_LED}(1:end-2)];

        signal_LED_time_parallel_output{i_LED} = reshape(signal_LED_time_reflection{i_LED}, [], Num_frame)';
        sig_ofdm_time{i_LED} = signal_LED_time_parallel_output{i_LED};
    end
end
