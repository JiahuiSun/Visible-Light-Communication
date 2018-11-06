function [estimated_channel_gain, estimated_distance] = curveFit(Num_LED, K, range, Hk_square)
detector_area = 1e-4; % unit: m^2
semiangle_at_half_power = pi/3;
parameter_m = -log(2)/log(cos(semiangle_at_half_power));
height_room = 5;
height_receiver = 1;

estimated_channel_gain = zeros(Num_LED, 1);
estimated_distance = zeros(Num_LED, 1);
% Curve fitting
myfit_type = fittype('A+B*cos(2*pi/128*x)+C*(2*cos(2*pi/128*x).^2-1)','independent','x','coefficients',{'A','B','C'});
myfit_opt = fitoptions(myfit_type);
myfit_opt.StartPoint = [ 0 0 0 ];  

for i_LED = 1:Num_LED
    K_part = [ K(range(i_LED,1):range(i_LED,2)), K(range(i_LED,3):range(i_LED,4)) ];
    Hk_square_part = [Hk_square(range(i_LED,1):range(i_LED,2)), Hk_square(range(i_LED,3):range(i_LED,4))];
    Hk_square_part_normalized = Hk_square_part./max(Hk_square_part);

    myfit = fit(K_part.', Hk_square_part_normalized.', myfit_type, myfit_opt);

%     figure();
%     plot(myfit,K_part,Hk_square_part_normalized);
%     hold on; grid on;
%     xlabel('K'); ylabel('|Hk|^2'); 
%     title(['Curve Fitting: ','LED ',num2str(i_LED)]);

    A = myfit.A;
    B = myfit.B;
    C = myfit.C;
    % prerequisite: h1>h3;
    D = ( sqrt(A+B+C) + sqrt(A-B+C) ) / 2;
    E = sqrt( D.^2 - 2*C );

    h1 = (D+E)/2;
    estimated_channel_gain(i_LED) = h1 * sqrt(max(Hk_square_part));
    estimated_distance(i_LED) = ( estimated_channel_gain(i_LED)  /((parameter_m + 1)*detector_area / (2*pi) * (height_room - height_receiver)^(parameter_m + 1))) ^ (-1/(parameter_m + 3));
end