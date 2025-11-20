% run_detection  主脚本：生成OFDM信号并执行盲检测算法
%
% 本脚本通过调用各个独立的函数，演示如何利用功率谱再处理方法
% 对未知OFDM信号的子载波间隔和子载波数量进行盲估计。
% 使用中文注释解释每一步的处理过程。

clear; close all; clc;

%% 参数设置
fs = 10e6;               % 采样频率10 MHz
num_subcarriers_true = 64;  % 真正的OFDM子载波数，用于生成信号
num_symbols = 10;          % 生成的OFDM符号个数
cp_ratio = 0.25;           % 循环前缀长度占有效符号长度的比例

%% 生成OFDM基带信号
[t, rx_signal, params] = generate_ofdm_signal(num_subcarriers_true, num_symbols, cp_ratio, fs);
fprintf('生成的OFDM信号长度: %d 样本\n', length(rx_signal));
fprintf('实际子载波数: %d, 理论子载波间隔: %.2f kHz\n', params.N, params.delta_f/1e3);

%% Step1: 利用Welch法估计功率谱
segment_length = 2048; % Welch窗口长度，可根据需要调整
[freq, Pxx] = compute_power_spectrum(rx_signal, fs, segment_length);

%% Step2: 估计信号带宽
threshold_ratio = 0.1; % 取功率谱最大值的10%作为阈值
[bandwidth, f_low, f_high] = estimate_bandwidth(freq, Pxx, threshold_ratio);
fprintf('估计信号带宽约: %.2f MHz (%.2f MHz ~ %.2f MHz)\n', bandwidth/1e6, f_low/1e6, f_high/1e6);

%% Step3: 对功率谱进行再处理
reprocessed = spectrum_reprocess(Pxx);

%% Step4: 估计再处理谱的周期
[M, peak_locs] = estimate_period(reprocessed);
fprintf('谱再处理序列主周期(点数): %d\n', M);

%% Step5: 根据带宽和周期估计子载波数
[N_est, delta_f_est] = estimate_subcarrier_number(bandwidth, fs, M);
fprintf('估计的子载波间隔: %.2f kHz\n', delta_f_est/1e3);
fprintf('估计的子载波数量: %d\n', N_est);

%% 可视化结果
figure;
subplot(3,1,1);
plot(freq/1e6, 10*log10(Pxx));
xlabel('频率 (MHz)'); ylabel('功率谱密度 (dB/Hz)');
title('平均功率谱'); grid on;
hold on;
% 标注估计的带宽区间
yl = ylim;
plot([f_low f_low]/1e6, yl, 'r--');
plot([f_high f_high]/1e6, yl, 'r--');
legend('功率谱','带宽下界','带宽上界');

subplot(3,1,2);
plot((0:length(reprocessed)-1), reprocessed);
xlabel('频率索引'); ylabel('再处理谱功率');
title('谱再处理结果'); grid on;
hold on;
% 标出检测到的峰值
plot(peak_locs, reprocessed(peak_locs), 'ro');
legend('谱再处理','检测到的峰');

subplot(3,1,3);
bar(1:length(peak_locs)-1, diff(peak_locs));
xlabel('峰序号'); ylabel('相邻峰间距(点)');
title(['估计周期M = ', num2str(M), '点']); grid on;

%% 打印总结
fprintf('\n===== 总结 =====\n');
fprintf('通过功率谱再处理方法估计OFDM信号参数:\n');
fprintf(' - 带宽约 %.2f MHz\n', bandwidth/1e6);
fprintf(' - 子载波间隔估计 %.2f kHz\n', delta_f_est/1e3);
fprintf(' - 子载波数估计为 %d，实际值 %d\n', N_est, num_subcarriers_true);
