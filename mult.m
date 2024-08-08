%%  信号参数设置
clear;clc
fs = 5e3;                    % 采样频率
fn = 1125;                    % 周期性冲击信号共振频率
fn2 = 2250;                   % 调制干扰信号共振频率
a = 100;                      % 衰减系数 
A0 = 0.08;                    % 周期性冲击信号位移常数
SNR10 = -10;                     %信噪比
SNR5 = -5;                     %信噪比
SNR3 = -3;                     %信噪比

N = 5*fs;                    % 采样点数
%%  调制干扰信号
tt=1/fs:1/fs:5;
s2 = 0.2*sin(pi*(20*tt));
sign_indices = find(sign(s2(1:end-1)) ~= sign(s2(2:end)));
index=sign_indices(1:2:end);
p1 = length(index);            % 重复次数
t=tt(index);
T = diff(t);
%%  周期性冲击信号
s12 = zeros(1,length(tt));
for i = 1:p1-1                               %产生p1个相同波形
    NT = round(fs*T(i));            % 单周期采样点数
    tt0 = 0:1/fs:(NT-1)/fs;      % 单周期采样时刻
    l=length(tt0);
    s12(index(i):index(i)+l-1) = (60*A0.*exp(-a.*(tt0)).*cos(2*pi*fn*(tt0))); 
end

s1 = find(s12 ~= 0);
s1 = s12(1:length(s1));
s2 = s2(1:length(s1));
%%  随机噪声
NOISE = randn(size(s1));
NOISE = NOISE-mean(NOISE);
signal_power = 1/length(s1)*sum(s1.*s1);
noise_variance10 = signal_power/(10^(SNR10/10));
noise_variance5 = signal_power/(10^(SNR5/10));
noise_variance3 = signal_power/(10^(SNR3/10));
NOISE10 = sqrt(noise_variance10)/std(NOISE)*NOISE;
NOISE5 = sqrt(noise_variance5)/std(NOISE)*NOISE;
NOISE3 = sqrt(noise_variance3)/std(NOISE)*NOISE;

%%  合成信号
s10 = s1+s2+NOISE10 ;
s5 = s1+s2+NOISE5 ;
s3 = s1+s2+NOISE3 ;

tt=1/fs:1/fs:length(s2)/fs;

%%  信号的时域图
figure(1)
plot(tt(1:length(s2)),s3)
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s)','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
xlim([0 3])
ylim([-5 5])
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);


wavename='cmor3-3';
totalscal=256;
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalscal;
scals=c./(1:totalscal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs=cwt(s1,scals,wavename); 

wlen = 80;
hop = wlen/8;
nfft = 4*wlen;

% generate analysis and synthesis windows
anal_win = blackmanharris(wlen, 'periodic');
synth_win = blackmanharris(wlen, 'periodic');



% perform time-frequency analysis and resynthesis of the signal
[STFT, f, t] = stft(s5, anal_win, hop, nfft, fs);
figure(5)
imagesc(t,f,abs(STFT));
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s) ','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Frequency(Hz)', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar


[u,a,v]=svd(STFT);
[m,n]=size(a);
A=zeros(m,n);
A(1,1)=a(1,1);
S=u*A*v';

singular_values = diag(a);

% 对奇异值进行排序
sorted_singular_values = sort(singular_values, 'descend');
% 提取前 50 个奇异值
top_50_singular_values = sorted_singular_values(1:20);

% 绘制柱状图
bar(top_50_singular_values);
xlabel('Index');
ylabel('Singular Value');
title('Top 20 Singular Values');
set(gca, 'FontName', 'Times New Roman');
[m,n]=size(STFT);


muzero=7e-2;    % the only tuning parameter

lambda=1e-5;  % model parameter
type=21;   %different modeling of Sparsity.
rate=0.05;   %update rate of \mu
gamma=0.01;     %gamma parameter in the rank approximation
tol=1e-4;  % stopping criterion

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
rk=rank(L1)

figure(6)
imagesc(t,f,abs(S));
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s) ','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Frequency(Hz)', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar



%%L1 norm
[M,N] = size(STFT);
lambda = 0.02 / sqrt(max(M,N)); 
u =  0.9*lambda;


% 调用 GPU 加速的 adm 函数
tic
[L2, S2] = adm(STFT, lambda, u, 500, 1e-1);
toc
rank(L2)
figure(9)
imagesc(t,f,abs(S2))
set(gcf,'Position',[20 100 660 460]);	 
xlabel('Time(s) ','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Frequency(Hz)', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
colorbar



[x_istft1, t_istft] = istft(L1, anal_win, synth_win, hop, nfft, fs);
[x_istft2, t_istft] = istft(L2, anal_win, synth_win, hop, nfft, fs);



%%VMD
alpha = fs/6;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
K = 3;              % 3 modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-5;

%--------------- Run actual VMD code

[u, u_hat, omega] = VMD(s10, alpha, tau, K, DC, init, tol);


%%FMD
filtersize = 30;
cutnum = 7;
modenum = 3;
maxiternum = 20;
% FMD
tic
y_final = FMD(16000, s10', filtersize, cutnum, modenum, maxiternum);
toc

figure(10)
plot(tt(1:length(s1)), s1, 'b', 'DisplayName', 'Noise-free signal')
hold on
plot(tt(1:length(x_istft2)), x_istft1, 'r', 'DisplayName', 'Our method')
plot(tt(1:length(x_istft2)), x_istft2, 'k', 'DisplayName', 'L1 method')
set(gcf,'Position',[20 100 1010 460]);	 
xlabel('Time(s)','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
xlim([0 3])
ylim([-5 7])
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
legend('Location', 'north', 'Interpreter', 'latex', 'Orientation', 'horizontal');

figure(11)
plot(tt(1:length(s1)), s1, 'b', 'DisplayName', 'Noise-free signal')
hold on
plot(tt(1:length(u)), u(2,:), 'g', 'DisplayName', 'VMD')
plot(tt(1:length(x_istft2)), x_istft1, 'r', 'DisplayName', 'Our method')
set(gcf,'Position',[20 100 1010 460]);	 
xlabel('Time(s)','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
xlim([0 3])
ylim([-5 7])
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
legend('Location', 'north', 'Interpreter', 'latex', 'Orientation', 'horizontal');

figure(12)
plot(tt(1:length(u)), y_final(:,2), 'g', 'DisplayName', 'FMD')
hold on
plot(tt(1:length(s1)), s1, 'b', 'DisplayName', 'Noise-free signal')
plot(tt(1:length(x_istft2)), x_istft1, 'r', 'DisplayName', 'Our method')
set(gcf,'Position',[20 100 980 460]);	 
xlabel('Time(s)','Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 12,'FontName','Times New Roman');
xlim([0 3])
ylim([-40 75])
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
legend('Location', 'north', 'Interpreter', 'latex', 'Orientation', 'horizontal');


x= x_istft2;
x_p = max(abs(x));  % 计算信号的峰值
mean_abs_x = mean(abs(x));  % 计算信号的绝对值的平均值
I = x_p / mean_abs_x;  % 计算脉冲指标

% 显示结果
disp(['信号的脉冲指标为：', num2str(I)]);


S = [4.8393, 5.9106, 7.4274, 8.6215;
     4.8364, 5.9533, 6.5643, 8.2295;
      5.0955, 6.2179, 6.6432, 6.6432];
x = 1:size(S, 1);  % 使用数字表示 x 轴位置
x_labels = {'-3dB', '-5dB', '-10dB'};

bar(x, S);  % 绘制条形图
ylim([0,9.8])
set(gca, 'XTickLabel', x_labels,'FontName','Times New Roman');  % 设置 x 轴标签
set(gcf,'Position',[20 100 980 460]);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', 24,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gcf,'Color','w');
set(gca, 'box', 'on', 'LineWidth', 2);
legend('VMD', 'FMD', 'L1 method','Ours method','FontSize', 20,'FontName','Times New Roman','Location', 'best','Orientation', 'horizontal'); % 添加图例并放置在最佳位置


kurtosis(u(2,:))
kurtosis(y_final(:,2))
kurtosis(x_istft1)
kurtosis(x_istft2)