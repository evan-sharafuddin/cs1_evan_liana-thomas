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

%% Processing the jazz signals
% loading audio -- space station has different sampling frequency for some
% reason
giant_steps = audioread("Giant Steps Bass Cut.wav");
[space_station,space_freq] = audioread("Space Station - Treble Cut.wav");
giant_time = 0:dt:(length(giant_steps)*dt-dt);
space_time = 0:dt:(length(space_station)*dt-dt);

% filtering audio -- need to filter both channels
giant_steps_bass = zeros(length(giant_steps),2);
space_station_treble = zeros(length(space_station),2);
giant_steps_bass(:,1) = bass(giant_steps(:,1),giant_time,dt);
giant_steps_bass(:,2) = bass(giant_steps(:,2),giant_time,dt);
space_station_treble(:,1) = treble(space_station(:,1),space_time,dt);
space_station_treble(:,2) = treble(space_station(:,2),space_time,dt);

% before and after
sound(giant_steps,freq),pause(3),clear sound
sound(giant_steps_bass,freq),pause(5),clear sound
sound(space_station,space_freq),pause(3),clear sound
sound(space_station_treble,space_freq),pause(5),clear sound

% visualizing the changes -- only looking at one channel
figure,hold on
subplot(2,1,1),spectrogram(giant_steps(:,1),1048,200,1048,freq)
title("Giant Steps, before bass-boosted preset"),clim([-140 -30]),xlim([0 12])
subplot(2,1,2),spectrogram(giant_steps_bass(:,1),1048,200,1048,freq)
title("Giant Steps, after bass-boosted preset"),clim([-140 -30]),xlim([0 12])
hold off

figure,hold on
subplot(2,1,1),spectrogram(space_station(:,1),1048,200,1048,space_freq);
title("Space Station, before treble boost preset"),clim([-140 -30]),xlim([0 12])
subplot(2,1,2),spectrogram(space_station_treble(:,1),1048,200,1048,space_freq)
title("Space Station, after treble boost preset"),clim([-140 -30]),xlim([0 12])
hold off

%%
figure,hold on
spectrogram(giant_steps(:,1),1048,200,1048,space_freq);
title("giant, before treble boost preset"),clim([-140 -30]),xlim([0 12])
figure,spectrogram(giant_steps_bass(:,1),1048,200,1048,space_freq)
title("gaint, after bass boost preset"),clim([-140 -30]),xlim([0 12])
hold off

[~,~,~,before_filter] = spectrogram(giant_steps(:,1),1048,200,1048,freq);
[~,t,f,after_filter] = spectrogram(giant_steps_bass(:,1),1048,200,1048,freq);

difference = abs(after_filter - before_filter);
figure,imagesc(t(end:-1:1), f, flipud( 20.*log10(difference) )' )
xlim([0 12e3])
colorbar


%%

% demonstrating the unity filter
giant_steps_unity = zeros(length(giant_steps),2);
giant_steps_unity(:,1) = unity(giant_steps(:,1),giant_time,dt);
giant_steps_unity(:,2) = unity(giant_steps(:,2),giant_time,dt);

figure,hold on
subplot(2,1,1),spectrogram(giant_steps(:,1),1048,200,1048,freq)
title("Giant Steps, before unity preset"),clim([-140 -30]),xlim([0 12])
subplot(2,1,2),spectrogram(giant_steps_unity(:,1),1048,200,1048,freq)
title("Giant Steps, after unity preset"),clim([-140 -30]),xlim([0 12])
hold off

%% removing unwanted background noise
blue_green = audioread("Blue in Green with Siren.wav");
blue_time = 0:dt:length(blue_green)*dt-dt;

% listen to before and process before
sound(blue_green,freq)

% siren noise occurs starting at ~ 1kHz, therefore, we want the first two
% bands of the equlaizer to be set to unity while the other three should be
% mostly/completely attenuated
blue_green_filtered = zeros(size(blue_green));
blue_green_filtered(:,1) = myFilter(blue_green(:,1),blue_time,dt,[2 2 0 0 0],5);
blue_green_filtered(:,2) = myFilter(blue_green(:,2),blue_time,dt,[2 2 0 0 0],5);

%%
% listening to sound after filter
sound(blue_green_filtered,freq)
figure,hold on
subplot(2,1,1),spectrogram(blue_green(:,1),1048,200,1048,freq)
title("Blue in Green with siren"),xlim([0 12])
subplot(2,1,2),spectrogram(blue_green_filtered(:,1),1048,200,1048,freq)
title("Blue in Green without siren"),xlim([0 12])
hold off

%% bird thing
bird_input = audioread("bird sound - Made with Clipchamp.mp4");
bird_time = 0:dt:length(bird_input)*dt-dt;

bird_filtered = zeros(size(bird_input));
bird_filtered(:,1) = myFilter(bird_input(:,1),bird_time,dt,[0 0 1 0 0],5);
bird_filtered(:,2) = myFilter(bird_input(:,2),bird_time,dt,[0 0 1 0 0],5);
figure,hold on
subplot(2,1,1),spectrogram(bird_input(:,1),1048,200,1048,freq),title("Bird recording before filter")
subplot(2,1,2),spectrogram(bird_filtered(:,1),1048,200,1048,freq),title("Bird recording after filter")
hold off

sound(bird_input,freq),pause(5),clear sound
sound(bird_filtered,freq)
% woo it works to isolate that one call pretty well at the end of the clip

%% LOAD AUDIO
filename = "Giant Steps Bass Cut.wav";
input = audioread(filename);
input = input(1:5/dt);

%% FILTER THAT AUDIOOOO

%FILTERRRRRRZZZZ

%gainzzzzz
% gain = [0 0 1 0 0]; % setting all gains equal to 1 seems to be pretty good unity setting
time = 0:dt:5-dt;
%can throw this into a loop later
% out = myFilter(input,time, dt,gain,5);

% coeff = [25, 200, 250, 500, 750, 1600, 2000, 15000];
% plot_all(coeff)

outTreble = treble(input, time, dt);
outBase = bass(input, time, dt);
outUnity = unity(input, time, dt);

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
sound(outTreble,freq),pause(5),clear sound
%%
% play sounds in each channel
% sound(out1,freq),pause(5),clear sound
% sound(out2,freq),pause(5),clear sound % at the moment, channels 1 and 2 are essentially the same 
% sound(out3,freq),pause(5),clear sound
% sound(out4,freq),pause(5),clear sound
% sound(out5,freq),pause(5),clear sound

%% FUNCTIONS BELOW

%% TREBLE BOOST
function output = treble (inputSig, time, dt)
gain = [1 1 1 5 5];
output = myFilter(inputSig,time, dt,gain,5);

end

%% BASS BOOST
function output = bass(audio_in, time, dt)
gain = [10 5 1 1 1];
output = myFilter(audio_in,time, dt,gain,5);
end

%% UNITY
function output = unity(audio_in, time, dt)
gain = [1 1 1 1 1]; % setting all gains equal to 1 seems to be pretty good unity setting
output = myFilter(audio_in,time, dt,gain,5);
end

%% Gain to dB
function out = gainToDb(gain_val)
out = 20*log(gain_val);
end


%% FILTER
function output = myFilter(input,time, dt,gain,iter)
if(length(gain) ~= 5)
    return
end
%DEFINIED CUTOFF FREQUENCIES
lp_fc = 150;
bp1_fc1 = 150;
bp1_fc2 = 300;
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
%% for the preset responses
function plot_preset ()
% inputs needs to be length 8
freq = 44.1e3;
dt = 1/freq;
time = 0:1/freq:4;

f= logspace(1.3, 4.3);

w = 2*pi*f;
response = zeros(size(f));

figure;
hold on;
% for every frequency in the logspace, create the input function and run it
% through a highpass and lowpass filter. Take the end time input to output
% ratio, as at the end the system should be stable.
for i = 1:3
    if (i == 1)
        for x=1:length(f)
            in = exp(1i*w(x)*time);
            Hout = bass(in, time, dt);
            response(x) = Hout(end) / in(end);
        end
        
    elseif(i == 2)
        for x=1:length(f)
            in = exp(1i*w(x)*time);
            Hout = unity(in, time, dt);
            response(x) = Hout(end) / in(end);
        end
    
    else
        for x=1:length(f)
            in = exp(1i*w(x)*time);
            Hout = treble(in, time, dt);
            response(x) = Hout(end) / in(end);
        end
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

legend("Bass", "Unity", "Treble");
end

%% Find the impulse responses
function impulse_response(coeff)
% inputs needs to be length 8
freq = 44.1e3;
dt = 1/freq;
time = 0:1/freq:.0075;

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
            
        in = zeros(size(time));
        in(1) = 1;
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
    
    subplot(2, 1, 1);
    hold on;
    plot(time, abs(Hout));
    title("|H(w)|");
    xlabel("Time (s)");
    ylabel("Magnitude");
    subplot(2, 1, 2);
    hold on;
    semilogx(time, angle(Hout));
    title("H(w) Angle");
    xlabel("Time");
    ylabel("Angle");

end
legend("Filter 1", "Filter 2", "Filter 3", "Filter 4", "Filter 5");

end

function impulse_preset ()
% inputs needs to be length 8
freq = 44.1e3;
dt = 1/freq;
time = 0:1/freq:.0075;

f= logspace(1.3, 4.3);

w = 2*pi*f;
response = zeros(size(f));

figure;
hold on;
% for every frequency in the logspace, create the input function and run it
% through a highpass and lowpass filter. Take the end time input to output
% ratio, as at the end the system should be stable.
for i = 1:3
    if (i == 1)
            in = zeros(size(time));
            in(1) = 1;
            Hout = bass(in, time, dt);
        
    elseif(i == 2)
            in = zeros(size(time));
            in(1) = 1;
            Hout = unity(in, time, dt);
    
    else
            in = zeros(size(time));
            in(1) = 1;
            Hout = treble(in, time, dt);
    end
    
    
    subplot(2, 1, 1);
    hold on;
    plot(time, abs(Hout));
    title("|H(w)|");
    xlabel("Time (s)");
    ylabel("Magnitude");
    subplot(2, 1, 2);
    hold on;
    semilogx(time, angle(Hout));
    title("H(w) Angle");
    xlabel("Time");
    ylabel("Angle");

end

legend("Bass", "Unity", "Treble");
end
