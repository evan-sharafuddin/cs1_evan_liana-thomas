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

%% FILTER THAT AUDIOOOO

%JUST TESTING each component filter
%out1 = lowpass(freq,input,200,dt);
%out1 = highpass(freq,input,8000,dt);
%out1 = bandpass(freq,input,100,500,5,dt);

%FILTERRRRRRZZZZ

%gainzzzzz
gain = [0 0 1 0 0]; % setting all gains equal to 1 seems to be pretty good unity setting

%can throw this into a loop later
out = myFilter(input,dt,gain,5);

% extracts the frequencies in each channel -- use to optomize the frequency
% cutoff values
out1 = myFilter(input,dt,[1 0 0 0 0],5);
out2 = myFilter(input,dt,[0 1 0 0 0],5);
out3 = myFilter(input,dt,[0 0 1 0 0],5);
out4 = myFilter(input,dt,[0 0 0 1 0],5);
out5 = myFilter(input,dt,[0 0 0 0 1],5);
%%
%TEST: Just play original audio
sound(input,freq),pause(5),clear sound

% PLAY filtered audios
sound(out,freq),pause(5),clear sound
%%

% play sounds in each channel
sound(out1,freq),pause(5),clear sound
sound(out2,freq),pause(5),clear sound % at the moment, channels 1 and 2 are essentially the same 
sound(out3,freq),pause(5),clear sound
sound(out4,freq),pause(5),clear sound
sound(out5,freq),pause(5),clear sound

%% FUNCTIONS BELOW

%% FILTER
function output = myFilter(input,dt,gain,iter)
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
y = gain(1) * 1.5*lowpass(input,lp_fc,dt,3); % LP filter iter controlled manually
%Bandpass
y = y + gain(2) * 1.5/atten_fact(bp1_fc1,bp1_fc2,iter,dt)*bandpass(input,bp1_fc1,bp1_fc2,iter,dt);
%Bandpass
y = y + gain(3) * 1.5/atten_fact(bp2_fc1,bp2_fc2,iter,dt)*bandpass(input,bp2_fc1,bp2_fc2,iter,dt);
%Bandpas
y = y + gain(4) * 1.5/atten_fact(bp3_fc1,bp3_fc2,iter,dt)*bandpass(input,bp3_fc1,bp3_fc2,iter,dt);
%Highpass
y = y + gain(5) * 1.5*highpass(input,hp_fc,dt,3); % HP filter iter controlled manually

output = y;
end



%% FILTER FUNCTIONS

%LOWPASS
function output = lowpass(x,fc,dt,iter)
t =0:dt:5-dt;
tau = 1/2/pi/fc;
y = x;
for i = 1:iter
y = lsim(tf(1/tau, [1, 1/tau]),y,t);
end
output = y;
end

%HIGHPASS
function output = highpass(x,fc,dt,iter)
t =0:dt:5-dt;
tau= 1/2/pi/fc;
% Filtered output using lsim
y = x;
for i = 1:iter
y = lsim(tf([1 0], [1, 1/tau]),y,t);
end
output = y;
end

%BANDPASS
function output = bandpass(x,lp_fc,hp_fc,iter,dt)
t =0:dt:5-dt;
y=x;
for n = 1:iter
    y = lowpass(y,lp_fc,dt,1);
    y = highpass(y,hp_fc,dt,1);
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
    input = exp(j*2*pi*f(i)*t);
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

