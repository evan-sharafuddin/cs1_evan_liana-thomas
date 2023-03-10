
%% SET UP
freq = 44.1e3;
dt = 1/freq;
f= logspace(1.3, 4.3);
w = 2*pi*f;
response = zeros(size(f));

%% LOAD AUDIO
filename = "Blue in Green with Siren.wav";
input = audioread(filename);
input = input(:, 1);
input = input(1:5/dt);

%% FILTER THAT AUDIOOOO

%JUST TESTING each component filter
%out1 = lowpass(freq,input,200,dt);
%out1 = highpass(freq,input,8000,dt);
%out1 = bandpass(freq,input,100,500,5,dt);

%FILTERRRRRRZZZZ

%gainzzzzz
gain = [0 0 0 0 10];
time = 0:dt:5-dt;

%can throw this into a loop later
out = myFilter(freq,input,time,gain);

coeff = [25, 200, 250, 500, 750, 1600, 2000, 15000];

plot_all(coeff);

%TEST: Just play original audio
%sound(input,freq)

%% PLAY filtered audios
sound(out,freq),pause(5),clear sound

%Plot Mag and Phase
% for m=1:length(f)
%     response(m) = out1(end) / input(end);
% end
%plot_freq(response)

%% FUNCTIONS BELOW

%% FILTER
function output = myFilter(fs,input,time,gain)
if(length(gain) ~= 5)
    return
end
%DEFINIED CUTOFF FREQUENCIES
lp_fc = 25;
bp1_fc1 = 200;
bp1_fc2 = 250;
bp2_fc1 = 500;
bp2_fc2 = 750;
bp3_fc1 = 1600;
bp3_fc2 = 2000;
hp_fc =15000;

%Lowpass
y = gain(1) .* lowpass(fs,input,lp_fc,time);
%Bandpass
y = y + (gain(2) .* bandpass(fs,input,bp1_fc1,bp1_fc2,5,time));
%Bandpass
y = y + (gain(3) .* bandpass(fs,input,bp2_fc1,bp2_fc2,5,time));
%Bandpass
y = y + (gain(4) .* bandpass(fs,input,bp3_fc1,bp3_fc2,5,time));
%Highpass
y = y + (gain(5) .* highpass(fs,input,hp_fc,time));

output = y;
end

%% FILTER FUNCTIONS

%LOWPASS
function output = lowpass(fs,x,fc,time)
% t =0:dt:5-dt;
tau = 1/2/pi/fc;
y = lsim(tf(1/tau, [1, 1/tau]),x,time);
output = y;
end

%HIGHPASS
function output = highpass(fs,x,fc,time)
% t =0:dt:5-dt;
tau= 1/2/pi/fc;
% Filtered output using lsim
y = lsim(tf([1 0], [1, 1/tau]),x,time);
output = y;
end

%BANDPASS
function output = bandpass(fs,x,lp_fc,hp_fc,iter,time)
% t =0:dt:5-dt;
for n = 1:iter
    intermediate = lowpass(fs,x,lp_fc,time);
    y = highpass(fs,intermediate,hp_fc,time);
end
output = y;

end

%% PLOTTING FUNCTIONS

%BODE PLOT
function plot_freq(y)
f= logspace(1.301, 4.301);
figure;
subplot(2, 1, 1);
semilogx(f, 20*log(abs(y)));
title("H(w) Mag");
xlabel("Frequency");
ylabel("Db");
subplot(2, 1, 2);
semilogx(f, angle(y)/pi);
title("H(w) Angle");
xlabel("Frequency");
ylabel("Angle");
end 

%Attenuation Factor
function output = atten_fact(input,filtered)

f = logspace(1,5,135);
H_mag = zeros(length(f),1);
H_phase = zeros(length(f),1);

for i = 1:length(f)
    H_mag(i) = abs(filtered(end)/input(end));
    H_phase(i) = angle(filtered(end)/input(end));
end

attenuation_factor = max(H_mag);

output = attenuation_factor;
end

function plot_all (coeff)
% inputs needs to be length 8
freq = 44.1e3;

time = 0:1/freq:4;

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
        tao1 = 1/2/pi/coeff(1); % 20
        tao2 = 0;
    elseif(i == 2)
        tao1 = 1/2/pi/coeff(2); % 150
        tao2 = 1/2/pi/coeff(3); % 250
    
    elseif(i==3)
        tao1 = 1/2/pi/coeff(4); % 500
        tao2 = 1/2/pi/coeff(5); % 750
    
    elseif(i==4)
        tao1 = 1/2/pi/coeff(6); % 1600
        tao2 = 1/2/pi/coeff(7); % 2000
        
    else
        % tao1 = 1/2/pi/19000;
        tao1 = 0;
        tao2 = 1/2/pi/coeff(8); % 20000
    end
            
    for x=1:length(f)
        in = exp(1i*w(x)*time);
        if (tao1 == 0)
            % Hout = lsim([1, 0],[1, 1/tao2], in, t); %highpass
            Hout = highpass(f, in, tao2,time);

        elseif(tao2 == 0)
            % Hout = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass
            Hout = lowpass(f, in, tao1, time);

        else
            % intermediate = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass

            % Hout = lsim([1, 0],[1, 1/tao2], intermediate, t); %highpass
            
            Hout = bandpass(f,in, tao1, tao2, 5,time);

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
end
