%% run.m  参数配置入口
% 将所有可调参数集中在此脚本，方便批量测试或记录实验配置。
% 运行本脚本即可完成配置并调用 main_ofdm_blind_det_ps。

cfg = struct();

%% 数据与限长配置
cfg.dataFile = 'starlink_ku_band_signal_480MHz.dat';       % 交错IQ文件 (int16 或 float32)
cfg.fs = 480e6;                    % 采样率 Hz
cfg.maxSamples = 0;            % 0 表示读取完整文件

%% Welch PSD 参数
cfg.welchLen = 8192;              % FFT 段长度
cfg.welchWindow = 'hann';         % 典型值: 'hann','hamming','rectwin'
cfg.welchOverlap = 0.5;           % 相邻段重叠比例 0~0.95

%% 再处理与检测参数
cfg.useImprovedReprocess = true;  % 启用平方根/瞬时幅度改进版
cfg.alphaIndex = 1;               % 改进版可用多个 α 结果, 1=默认
cfg.zeroPadFactor = 2;            % 基本版功率谱再处理的零填充倍率
cfg.suppressZeroLagPeak = true;   % 是否抑制再处理结果的零延迟巨大峰值
cfg.periodicityThreshold = 2.0;   % 周期性指标阈值
cfg.minPeakProminence = 0.02;     % 周期检测最小峰值显著度
cfg.smoothingSpan = 51;           % 周期检测包络平滑窗口长度
cfg.bandwidthThreshold_dB = 6;    % 带宽检测相对阈值 (dB)
cfg.bandwidthGuardFraction = 0.05;% 噪声底估计用的两侧保护带比例
cfg.bandwidthEnergyCoverage = 0.995; % 采用能量百分位缩小带宽 (0~1)
cfg.bandwidthDetectionMode = 'hybrid'; % 'threshold'|'energy'|'hybrid'
cfg.bandwidthMinWidthHz = [];     % 覆盖最小带宽 (Hz)，为空则自动=Δf

%% 可视化与显示控制
cfg.plotFigures = true;           % 生成所有图形
cfg.closeAllFigures = true;       % 运行前关闭所有图形窗口
cfg.clearCommandWindow = true;    % 运行前清理命令行

%% 调用主流程
results = main_ofdm_blind_det_ps(cfg);

disp('--- 运行完成，可查看 results 结构体了解详细输出 ---');
