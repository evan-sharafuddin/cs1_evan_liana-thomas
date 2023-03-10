t = 0:1/freq:4;

tao1 = 0;
tao2 = 0;


f= logspace(1.3, 4.3);

w = 2*pi*f;

response = zeros(size(f));

figure;
hold on;
% for every frequency in the logspace, create the input function and run it
% through a highpass and lowpass filter. Take the end time input to output
% ratio, as at the end the system should be stable.
for i = 1:5
    if (i == 1)
        tao1 = 1/2/pi/20;
        tao2 = 0;
    elseif(i == 2)
        tao1 = 1/2/pi/150;
        tao2 = 1/2/pi/200;
    
    elseif(i==3)
        tao1 = 1/2/pi/500;
        tao2 = 1/2/pi/750;
    
    elseif(i==4)
        tao1 = 1/2/pi/1600;
        tao2 = 1/2/pi/2000;
        
    else
        % tao1 = 1/2/pi/19000;
        tao1 = 0;
        tao2 = 1/2/pi/20000;
    end
            
    for x=1:length(f)
        in = exp(1i*w(x)*t);
        if (tao1 == 0)
            Hout = lsim([1, 0],[1, 1/tao2], in, t); %highpass

        elseif(tao2 == 0)
            Hout = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass

        else
            intermediate = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass

            Hout = lsim([1, 0],[1, 1/tao2], intermediate, t); %highpass

        end
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