% I got the filters to filter a lot better. The cutoff frequencies still need to be optomized 
% though and it might be nice for Thomas to graph some of the frequency responses of some 
% of the filters. 

% in myFilters() I multiplied the output from each filter by values to get them back 
% up to the original audio level. I multipled the lowpass and highpass by 1.5 (trial
% and error) but they still seemed quiet so might be good to play around with that. 
% For the bandpass filters I multiplied by the reciprocal of the attenuation factor 
% thing to kind of "cancel out" the attenuation, but then also multiplied by 1.5 to 
% get it a little bit louder. Might be nice to change this as you change the 
% cutoff frequencies but with everything set up the unity preset works pretty well

% Let me know if anything doesn't make sense --Evan




%% SET UP
freq = 44.1e3;
dt = 1/freq;
f= logspace(1.3, 4.3);
w = 2*pi*f;
response = zeros(size(f));

%% LOAD AUDIO
filename = "Blue in Green with Siren.wav";
input = audioread(filename);

input = input(1:5/dt);

%% Filtering Out Radio Static
am_input = audioread("bird_sound - Made with Clipchamp.mp4");
am_input = am_input(1:5/dt);

figure, spectrogram(am_input,1048,200,1048,freq),title("Radio before filter")

am_filtered = myFilter(am_input,0:dt:5-dt, dt,[0 0 1 0 0],5);

figure, spectrogram(am_filtered,256,200,256,freq),title("Radio after filter")

% woo it works to isolate that one call pretty well at the end of the clip

%% FILTER THAT AUDIOOOO

%FILTERRRRRRZZZZ

%gainzzzzz
% gain = [0 0 1 0 0]; % setting all gains equal to 1 seems to be pretty good unity setting
time = 0:dt:5-dt;
%can throw this into a loop later
% out = myFilter(input,time, dt,gain,5);

% coeff = [25, 200, 250, 500, 750, 1600, 2000, 15000];
% plot_all(coeff)

out = treble(input, time, dt);

% extracts the frequencies in each channel -- use to optomize the frequency
% cutoff values
% out1 = myFilter(input,time, dt,[1 0 0 0 0],5);
% out2 = myFilter(input,time,dt,[0 1 0 0 0],5);
% out3 = myFilter(input,time,dt,[0 0 1 0 0],5);
% out4 = myFilter(input,time,dt,[0 0 0 1 0],5);
% out5 = myFilter(input,time,dt,[0 0 0 0 1],5);
%%
%TEST: Just play original audio
sound(input,freq),pause(5),clear sound

% PLAY filtered audios
sound(out,freq),pause(5),clear sound
%%
% play sounds in each channel
% sound(out1,freq),pause(5),clear sound
% sound(out2,freq),pause(5),clear sound % at the moment, channels 1 and 2 are essentially the same 
% sound(out3,freq),pause(5),clear sound
% sound(out4,freq),pause(5),clear sound
% sound(out5,freq),pause(5),clear sound

%% FUNCTIONS BELOW
function output = treble (inputSig, time, dt)

output = myFilter(inputSig,time, dt,[0.2 0.4 1 2 2],5);

end

%% FILTER
function output = myFilter(input,time, dt,gain,iter)
if(length(gain) ~= 5)
    return
end
%DEFINIED CUTOFF FREQUENCIES
lp_fc = 150;
bp1_fc1 = 150;
bp1_fc2 = 200;
bp2_fc1 = 500;
bp2_fc2 = 1000;
bp3_fc1 = 3000;
bp3_fc2 = 5000;
hp_fc = 5000;

% Multiplying each of the filter outputs by 1.5 as well as the reciprocal
% of the attenuation factor for the bandpass filters seemed to work pretty
% well in restoring the original sound for the unity setting (i.e. leaving
% all of the gain values to be 1 and passing the audio through the filter
% resulted in mostly no change)

%Lowpass
y = gain(1) * 1.5*lowpass(input,lp_fc,time,3); % LP filter iter controlled manually
%Bandpass
y = y + gain(2) * 1.5/atten_fact(bp1_fc1,bp1_fc2,iter,dt)*bandpass(input,bp1_fc1,bp1_fc2,iter,time);
%Bandpass
y = y + gain(3) * 1.5/atten_fact(bp2_fc1,bp2_fc2,iter,dt)*bandpass(input,bp2_fc1,bp2_fc2,iter,time);
%Bandpas
y = y + gain(4) * 1.5/atten_fact(bp3_fc1,bp3_fc2,iter,dt)*bandpass(input,bp3_fc1,bp3_fc2,iter,time);
%Highpass
y = y + gain(5) * 1.5*highpass(input,hp_fc,time,3); % HP filter iter controlled manually

output = y;
end



%% FILTER FUNCTIONS

%LOWPASS
function output = lowpass(x,fc,time,iter)
tau = 1/2/pi/fc;
y = x;
for i = 1:iter
y = lsim(tf(1/tau, [1, 1/tau]),y,time);
end
output = y;
end

%HIGHPASS
function output = highpass(x,fc,time,iter)
tau= 1/2/pi/fc;
% Filtered output using lsim
y = x;
for i = 1:iter
y = lsim(tf([1 0], [1, 1/tau]),y,time);
end
output = y;
end

%BANDPASS
function output = bandpass(x,lp_fc,hp_fc,iter,time)
y=x;
for n = 1:iter
    y = lowpass(y,lp_fc,time,1);
    y = highpass(y,hp_fc,time,1);
end
output = y;

end

%% PLOTTING FUNCTIONS

%BODE PLOT
function plot_freq()
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

%Attenuation Factor -- basically does all the work of creating the plot --
%calculates the magnitude of the frequency response. The maximum magnitude
%is found and set as the "attenuation factor" (in other words, the selected
%frequency is attenuated by a factor of ___). By multiplying the output of
%hte bandpass filter with the reciprocal of this value, we are able to
%mostly get back to the original audio amplitude.

function output = atten_fact(lp_fc,hp_fc,iter,dt)
f = logspace(1,5,135);
H_mag = zeros(length(f),1);

lp_tau = 1/2/pi/lp_fc;
hp_tau = 1/2/pi/hp_fc;

% arbitrary time range but long enough to ensure steady state has been
% reached
t = 0:dt:100*hp_tau;

for i = 1:length(f)
    input = exp(1i*2*pi*f(i)*t);
    output = input;
 
    % simulating the bandpass filter here
    for x = 1:iter
        output = lsim(tf(1/lp_tau, [1, 1/lp_tau]),output,t);
        output = lsim(tf([1 0], [1, 1/hp_tau]),output,t);
    end
    H_mag(i) = abs(output(end));
end

attenuation_factor = max(H_mag);
output = attenuation_factor;
end

function plot_all (coeff)
% inputs needs to be length 8
freq = 44.1e3;
dt = 1/freq;
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
        tao1 = coeff(1); % 20
        tao2 = 0;
    elseif(i == 2)
        tao1 = coeff(2); % 150
        tao2 = coeff(3); % 250
    
    elseif(i==3)
        tao1 = coeff(4); % 500
        tao2 = coeff(5); % 750
    
    elseif(i==4)
        tao1 = coeff(6); % 1600
        tao2 = coeff(7); % 2000
        
    else
        % tao1 = 1/2/pi/19000;
        tao1 = 0;
        tao2 = coeff(8); % 20000
    end
            
    for x=1:length(f)
        in = exp(1i*w(x)*time);
        if (tao1 == 0)
            % Hout = lsim([1, 0],[1, 1/tao2], in, t); %highpass
            Hout = 1.5*highpass(in, tao2, time, 5);

        elseif(tao2 == 0)
            % Hout = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass
            Hout = 1.5*lowpass(in, tao1, time, 5);

        else
            % intermediate = lsim(1/tao1, [1, 1/tao1], in, t); % lowpass

            % Hout = lsim([1, 0],[1, 1/tao2], intermediate, t); %highpass
            
            Hout = 1.5/atten_fact(tao1,tao2,5,dt) * bandpass(in, tao1, tao2, 5, time);

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