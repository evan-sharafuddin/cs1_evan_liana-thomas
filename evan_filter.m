%% Filter attempt
clear, close all

fs = 44.1e3; % increased because otherwise systems exhibit DT behavior
dt = 1/fs;

% parameters
target_low = 500;
target_high = 550;
iter = 7;
% change iterations for lowpass and highpass filters
lp_iter = 10; hp_iter = 10;
gain = 1;

% defining circuit components
lp_Tao = 1/(2*pi*target_low);
hp_Tao = 1/(2*pi*target_high);
t = 0:dt:100*lp_Tao;


% visual of what the filter looks like 
figure
attenuation_factor = generate_bode(lp_Tao,hp_Tao,t,iter,0,0) % bandpass (LP & HP)
generate_bode(lp_Tao,0,t,lp_iter,1,0); % LP
generate_bode(0,hp_Tao,t,hp_iter,0,1); % HP
hold on 
legend('Bandpass','Lowpass','Highpass')
xlim([0 2*target_high])
%%

dt = 1/44.1e3;
input = audioread("rickroll.wav");
input = input(1:5/dt);
sound(input,44.1e3),pause(5),clear sound
output = bandpass(input,0:dt:5-dt,lp_Tao,hp_Tao,iter);

% amplify
output = output.*1/attenuation_factor*gain;
sound(output,44.1e3),pause(5),clear sound

%%
function output = bandpass(input,t,lp_Tao,hp_Tao,iter)
% runs lowpass and highpass systems in serial for the desired amount of
% iterations
output = input;
for x = 1:iter
    output = lowpass(output,t,lp_Tao);
    output = highpass(output,t,hp_Tao);
end  
end
%%
function output = lowpass(input,t,Tao)
    output = lsim(tf(1/Tao, [1, 1/Tao]),input,t);
end

function output = highpass(input,t,Tao)
    output = lsim(tf([1 0], [1, 1/Tao]),input,t);
end

%%
function attenuation_factor = generate_bode(lp_Tao,hp_Tao,t,iter,lp,hp)
f = logspace(1,5,135);
H_mag = zeros(length(f),1);
H_phase = zeros(length(f),1);
for i = 1:length(f)
    w = 2*pi*f(i);
    input = exp(j*w*t);

    output = input;
    
    for x = 1:iter
        % running just lp or hp
        if lp
            output = lowpass(output,t,lp_Tao);
        elseif hp
            output = highpass(output,t,hp_Tao);
            
        % run desired systems in serial
        else
            output = lowpass(output,t,lp_Tao);
            output = highpass(output,t,hp_Tao);
        end
    end

    % gathers frequency response for that particular frequency in magnitude
    % and phase
    H_mag(i) = abs(output(end)/input(end));
    H_phase(i) = angle(output(end)/input(end));
end

attenuation_factor = max(H_mag);

% graphing stuff
hold on
plot(f,H_mag)
title('|H(\omega)|')
xlabel('\omega/2\pi (Hz)')
ylabel('|H(\omega)|')
grid on
hold off

end
