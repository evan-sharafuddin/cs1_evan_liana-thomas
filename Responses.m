% humans can hear from 20 Hz to 20 kHz
% freq = 44.1e3;
% t = 0:1/freq:4;

% tao1 is for lowpass, higher value
% tao2 is for highpass, lower value

% cutoff should be 1/2pi/tao

% near 0
% tao1 = 1/2/pi;
% tao2 = 1/2/pi/100;

tao1 = 1/2/pi/500;
tao2 = 1/2/pi/3000;

% near 500
% tao1 = .55e-3;
% tao2 = .2e-3; % increasing this increases entire amplitude, moves peak back

% 1000 freq center
% tao1 = 1.6e3 * .1e-6;
% tao2 = 680*.2e-6;

% 5000
% tao1 = 1/2/pi/3000;
% tao2 = 1/2/pi/10000;

% near 20k
% tao1 = 1/2/pi/10000;
% tao2 = 1/2/pi/20000;

f= logspace(1.3, 4.3);

w = 2*pi*f;
% log of 20 is close to 1.3, log of 20k is close to 4.3
% responselow = zeros(size(f));
response = zeros(size(f));

figure;
hold on;
% for every frequency in the logspace, create the input function and run it
% through a highpass and lowpass filter. Take the end time input to output
% ratio, as at the end the system should be stable.
for i = 1:5
    if (i == 1)
        tao1 = 1/2/pi;
        tao2 = 1/2/pi/100;
    elseif(i == 2)
        tao1 = 1/2/pi/100;
        tao2 = 1/2/pi/500;
    
    elseif(i==3)
        tao1 = 1/2/pi/500;
        tao2 = 1/2/pi/3000;
    
    elseif(i==4)
        tao1 = 1/2/pi/3000;
        tao2 = 1/2/pi/10000;
        
    else
        tao1 = 1/2/pi/10000;
        tao2 = 1/2/pi/20000;
    end
            
    for x=1:length(f)
        in = exp(1i*w(x)*t);
        intermediate = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass
        Hout = lsim([1, 0],[1, 1/tao2], intermediate, t); %highpass
        response(x) = Hout(end) / in(end);

    end
    
    subplot(2, 1, 1);
    hold on;
    set(gca,'XScale','log')
    semilogx(f, 20*log(abs(response)));
    title("|H(w)|");
    xlabel("Frequency");
    ylabel("Db");
    subplot(2, 1, 2);
    hold on;
    set(gca,'XScale','log')
    semilogx(f, angle(response)/pi);
    title("H(w) Angle");
    xlabel("Frequency");
    ylabel("Angle");

end

legend("Filter 1", "Filter 2", "Filter 3", "Filter 4", "Filter 5");