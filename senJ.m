% %%  信号参数设置
clear;clc
% fs = 5e3;                    % 采样频率
% fn = 1125;                    % 周期性冲击信号共振频率
% fn2 = 2250;                   % 调制干扰信号共振频率
% a = 100;                      % 衰减系数 
% A0 = 0.08;                    % 周期性冲击信号位移常数
% SNR = -10;                     %信噪比
% 
% N = 5*fs;                    % 采样点数
% %%  调制干扰信号
% tt=1/fs:1/fs:5;
% s3 = 0.2*sin(pi*(20*tt));
% sign_indices = find(sign(s3(1:end-1)) ~= sign(s3(2:end)));
% index=sign_indices(1:2:end);
% p1 = length(index);            % 重复次数
% t=tt(index);
% T = diff(t);
% %%  周期性冲击信号
% s12 = zeros(1,length(tt));
% for i = 1:p1-1                               %产生p1个相同波形
%     NT = round(fs*T(i));            % 单周期采样点数
%     tt0 = 0:1/fs:(NT-1)/fs;      % 单周期采样时刻
%     l=length(tt0);
%     s12(index(i):index(i)+l-1) = (60*A0.*exp(-a.*(tt0)).*cos(2*pi*fn*(tt0))); 
% end
% 
% s1 = find(s12 ~= 0);
% s1 = s12(1:length(s1));
% s3 = s3(1:length(s1));
% %%  随机噪声
% NOISE = randn(size(s1));
% NOISE = NOISE-mean(NOISE);
% signal_power = 1/length(s1)*sum(s1.*s1);
% noise_variance = signal_power/(10^(SNR/10));
% NOISE = sqrt(noise_variance)/std(NOISE)*NOISE;
% 
% %%  合成信号
% s = s1+s3+NOISE ;

a=load('data.mat');
s=a.s;
s1=a.s1;
s3=a.s3;
NOISE=a.NOISE;
fs=a.fs;
tt=1/fs:1/fs:length(s3)/fs;

%%  信号的时域图
figure(1)
plot(tt(1:length(s3)),s)
hold on
plot(tt(1:length(s3)),s1,'r')
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s)','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
xlim([0 3])
ylim([-11 11])
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);


f0 = .01*fs;            % cycle frequency (in Hz)
Nw = 2^7;               % window length (number of samples)
alpha_max = 8*f0;       % maximum cyclic frequency to scan (in Hz)
                        % (should cover the cyclic frequency range of interest)
opt.coh = 1;            % compute sepctral coherence? (yes=1, no=0)
[CS,alpha,~,~] = Fast_SC(s,Nw,alpha_max,fs,opt);


figure(1)
% 绘制第一个子图
subplot(2, 1, 1)
plot(tt(1:length(s3)), s)
hold on
plot(tt(1:length(s3)), s1, 'r')
xlabel('Time(s)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman')
xlim([0 3])
ylim([-11 11])
set(gca, 'YDir', 'normal', 'FontSize', 24, 'box', 'on', 'LineWidth', 2)


% 绘制第二个子图
subplot(2, 1, 2)
plot(alpha(2:end), mean(abs(CS(:, 2:end))),'DisplayName', 'CS')
hold on
for x = 10:10:200
    plot([x x], [0.02 0.08], 'r--')
end
plot([x x], [0.02 0.08], 'r--', 'DisplayName', 'H')
xlabel('Cyclic frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times New Roman')
xlim([0 200])
ylim([0.03 0.08])
set(gca, 'YDir', 'normal', 'FontSize', 24, 'box', 'on', 'LineWidth', 2)

% 设置图形整体属性
set(gcf, 'Position', [20 100 660 660], 'Color', 'w')



wlen = 80;
hop = wlen/8;
nfft = 4*wlen;

% generate analysis and synthesis windows
anal_win = blackmanharris(wlen, 'periodic');
synth_win = blackmanharris(wlen, 'periodic');

% perform time-frequency analysis and resynthesis of the signal
[STFT, f, t] = stft(s, anal_win, hop, nfft, fs);
figure(2)
imagesc(t,f,abs(STFT));
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s) ','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Frequency(Hz)', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar 

[m,n]=size(STFT);


muzero=1.1;    % the only tuning parameter

lambda=4.5e-6;  % model parameter
type=21;   %different modeling of Sparsity.
rate=0.5;   %update rate of \mu
gamma=4e-4;     %gamma parameter in the rank approximation
tol=1e-3;  % stopping criterion

%initializations
S=zeros(m,n);
Y=zeros(m,n);
L1=STFT;
si=zeros(min(m,n),1); % for DC 
mu=muzero;

tic;
for ii=1:500
  
    D=STFT-S-Y/mu;
    [L1,si] = DC(D,mu/2,si,gamma);
    [S]=errorsol(Y,STFT,L1,lambda,mu,type);
    Y=Y+mu*(L1-STFT+S);
    mu=mu*rate;
    
    sigma=norm(STFT-S-L1,'fro');
    RRE=sigma/norm(STFT,'fro');
    
    if RRE<tol
        break
    end
    
end
time_cost = toc;


%%L1 norm
[M,N] = size(STFT);
lambda = 0.02 / sqrt(max(M,N)); 
u =  0.9*lambda;


% 调用 GPU 加速的 adm 函数
tic
[L2, S2] = adm(STFT, lambda, u, 500, 1e-1);
toc
rank(L2)


figure(3);
% 第一张子图（第一行居中）
subplot(3,2,[1,2]);
imagesc(t,f,abs(STFT)); 
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar

% 第二张子图（第二行左侧）
subplot(3,2,3);
imagesc(t,f,abs(L2)) 
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar
% 在此添加第二张子图的绘图代码

% 第三张子图（第二行右侧）
subplot(3,2,4);
imagesc(t,f,abs(L1))	 
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar

% 第四张子图（第三行左侧）
subplot(3,2,5);
imagesc(t,f,abs(S2))	 
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar

% 第五张子图（第三行右侧）
subplot(3,2,6);
imagesc(t,f,abs(S))	 
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar
% 在此添加第五张子图的绘图代码

set(gcf, 'Color', 'w');
left = 100;      % 左边距
bottom = 100;    % 底边距
width = 1200;    % 图片宽度
height = 1300;    % 图片高度
set(gcf, 'Position', [left bottom width height]);
annotation('textbox', [0.45, 0.03, 0.05, 0.05], 'String', 'Time(s)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');
annotation('textbox', [0.07, 0.4, 0.05, 0.05], 'String', 'Frequency(Hz)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none','Rotation', 90);
annotation('textbox', [0.3, 0.63, 0.05, 0.05], 'String', '=', 'Interpreter', 'latex', 'FontSize', 48, 'FontName', 'Times New Roman', 'EdgeColor', 'none','Rotation', 90);
annotation('textbox', [0.73, 0.63, 0.05, 0.05], 'String', '=', 'Interpreter', 'latex', 'FontSize', 48, 'FontName', 'Times New Roman', 'EdgeColor', 'none','Rotation', 90);
annotation('textbox', [0.3, 0.34, 0.05, 0.05], 'String', '+', 'Interpreter', 'latex', 'FontSize', 48, 'FontName', 'Times New Roman', 'EdgeColor', 'none','Rotation', 90);
annotation('textbox', [0.73, 0.34, 0.05, 0.05], 'String', '+', 'Interpreter', 'latex', 'FontSize', 48, 'FontName', 'Times New Roman', 'EdgeColor', 'none','Rotation', 90);
annotation('ellipse',[0.135 0.185 0.28 0.04], 'Color', 'red','LineWidth',2) 
annotation('ellipse',[0.575 0.185 0.28 0.04], 'Color', 'red','LineWidth',2) 
annotation('ellipse',[0.13 0.78 0.73 0.05], 'Color', 'red','LineWidth',2) 
x_line_between_1_and_2 = [0.49, 0.49];
y_line_between_1_and_2 = [0.1, 0.66];
annotation('line', x_line_between_1_and_2, y_line_between_1_and_2, 'Color', 'black', 'LineWidth', 3);
annotation('textbox', [0.074, 0.88, 0.05, 0.05], 'String', '(a)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');
annotation('textbox', [0.074, 0.59, 0.05, 0.05], 'String', '(b)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');
annotation('textbox', [0.074, 0.29, 0.05, 0.05], 'String', '(c)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');
annotation('textbox', [0.5, 0.59, 0.05, 0.05], 'String', '(d)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');
annotation('textbox', [0.5, 0.29, 0.05, 0.05], 'String', '(e)', 'Interpreter', 'latex', 'FontSize', 36, 'FontName', 'Times New Roman', 'EdgeColor', 'none');






