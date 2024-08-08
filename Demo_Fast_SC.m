% This script demonstrates the use of the Fast Spectral Correlation code
% "Fast_SC.m" on a synthetic cyclostationary signal. 

clear all
close all

% Synthesis of an elementary cyclostationary signal
% =================================================
Fs = 1e3;               % sampling frequency (in Hz)
L = 1e5;                % signal length (number of samples)
f0 = .01*Fs;            % cycle frequency (in Hz)

x = randn(L,1);

% filtration by a resonance
a = [1 -2*cos(2*pi*.2)*.9 .9^2];
x = filter(1,a,x);      

% periodic amplitude modulation
x = x.*(1 + sin(2*pi*(0:L-1)'*f0/Fs));

% addition of white noise (SNR = 0dB)
x = x + std(x)*randn(L,1);  

figure
plot((0:L-1)/Fs,x),title('Synthetic cyclostationary signal')
xlabel('time (s)')


% Fast Spectral Correlation
% =========================
% Analysis parameters to be fixed by the user
Nw = 2^7;               % window length (number of samples)
alpha_max = 4*f0;       % maximum cyclic frequency to scan (in Hz)
                        % (should cover the cyclic frequency range of interest)
opt.coh = 1;            % compute sepctral coherence? (yes=1, no=0)


[S,alpha,f,Nv] = Fast_SC(s,Nw,alpha_max,Fs,opt);


figure,subplot(211)
imagesc(alpha(2:end),f,abs(S(:,2:end))),axis xy,colorbar('Location','north')
if opt.coh == 0,title('Spectral Correlation'),else,title('Spectral Coherence'),end
xlabel('cyclic frequency \alpha (Hz)'),ylabel('f (Hz)'),
xlim([0 alpha_max])

subplot(212),plot(alpha(2:end),mean(abs(S(:,2:end))),'k'),xlim([0 alpha_max])
title('Enhanced envelope spectrum')
xlabel('cyclic frequency \alpha (Hz)')


% It may happen that the resolution of the computer screen is not enough to
% display all the details in the spectral correlation; in this case, a solution
% is to decrease the cycli spectral resolution by smoothing.

screen_resolution = 1080;       % number of vertical lines
Na = length(alpha);
N = round(Na/screen_resolution);
w = hamming(2*N+1);
if N > 1
    Smooth = zeros(Nw/2+1,Na);
    for k = 1:Nw/2+1
        Smooth(k,:) = conv(abs([0 S(k,2:end)]),w,'same');
    end
    
    figure,
    imagesc(alpha(2:end),f,Smooth),axis xy,colorbar('Location','north')
    if opt.coh == 0,title('Smoothed spectral Correlation'),else,title('Smoothed spectral Coherence'),end
    xlabel('cyclic frequency \alpha (Hz)'),ylabel('f (Hz)'),
    xlim([0 alpha_max])
end

plot(alpha)


