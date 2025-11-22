function results = main_ofdm_blind_det_ps(cfg)
%% main_ofdm_blind_det_ps  使用功率谱再处理法分析真实IQ数据的子载波参数
% 使用说明:
%   1) 将 cfg.dataFile 指向实际IQ文件 (交错IQ, int16或float32)。
%   2) 设置与采样系统一致的 cfg.fs (Hz)。
%   3) 调用本函数或运行 run.m，自动完成读取→功率谱→再处理→周期检测→参数估计。
%      所有关键步骤均以函数形式实现，可单独调用。

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
cfg = apply_main_defaults(cfg);

if cfg.closeAllFigures
    close all;
end
if cfg.clearCommandWindow
    clc;
end

dataFile = cfg.dataFile;
fs = cfg.fs;
maxSamples = cfg.maxSamples;
welchLen = cfg.welchLen;
welchWindow = cfg.welchWindow;
welchOverlap = cfg.welchOverlap;
useImprovedReprocess = cfg.useImprovedReprocess;
alphaIndex = cfg.alphaIndex;
periodicityThreshold = cfg.periodicityThreshold;

fprintf('=== 基于功率谱再处理的OFDM盲检测 (文件: %s) ===\n', dataFile);

%% 1. 读取IQ文件
[rx_signal_full, iqMeta] = read_iq_autodetect(dataFile, struct('format','auto','Verbose',true));
fprintf('检测到格式: %s, endian=%s, 顺序=%s\n', iqMeta.detected_format, iqMeta.endian, iqMeta.order);

if isempty(rx_signal_full)
    error('无法从 %s 读取到有效IQ数据。', dataFile);
end

if maxSamples > 0 && numel(rx_signal_full) > maxSamples
    fprintf('信号长度 %d 点, 仅截取前 %d 点参与分析。\n', numel(rx_signal_full), maxSamples);
    rx_signal = rx_signal_full(1:maxSamples);
else
    rx_signal = rx_signal_full;
end

% 去直流，便于功率谱估计
rx_signal = rx_signal - mean(rx_signal);

%% 2. Welch功率谱
[Pxx, f_grid, welchInfo] = estimate_psd_welch(rx_signal, fs, welchLen, ...
    'Window', welchWindow, 'Overlap', welchOverlap, 'CenterDC', true);
df_psd = welchInfo.delta_f;

%% 3. 带宽估计
bandwidthArgs = {'Threshold_dB', cfg.bandwidthThreshold_dB, ...
    'GuardFraction', cfg.bandwidthGuardFraction, ...
    'EnergyCoverage', cfg.bandwidthEnergyCoverage, ...
    'DetectionMode', cfg.bandwidthDetectionMode, ...
    'FigureTitle', dataFile + " - 带宽估计", 'PlotFigure', cfg.plotFigures};
if ~isempty(cfg.bandwidthMinWidthHz)
    bandwidthArgs = [bandwidthArgs, {'MinWidthHz', cfg.bandwidthMinWidthHz}];
end
[B_hat, f_L, f_H, ~] = estimate_bandwidth_from_psd(Pxx, f_grid, bandwidthArgs{:});

%% 4. 功率谱再处理 (基本版 + 可选改进版)
[R_basic, axis_basic] = spectrum_reprocessing_basic(Pxx, 'ZeroPadFactor', cfg.zeroPadFactor, ...
    'SuppressZeroLag', cfg.suppressZeroLagPeak, ...
    'FigureTitle', dataFile + " - 基本功率谱再处理", 'PlotFigure', cfg.plotFigures);

R_use = R_basic;
if useImprovedReprocess
    R_imp = spectrum_reprocessing_improved(Pxx, 'TimeSignal', rx_signal, 'Fs', fs, ...
        'FigureTitle', dataFile + " - 改进功率谱再处理", 'PlotFigure', cfg.plotFigures);
    if ~isempty(R_imp.root)
        alphaIdx = min(alphaIndex, numel(R_imp.root));
        R_use = R_imp.root{alphaIdx};
        reprocessDesc = sprintf('改进版 α=%.2f', R_imp.alphas(alphaIdx));
    else
        reprocessDesc = '基本版 (改进失败)';
    end
else
    reprocessDesc = '基本版';
end

%% 5. 周期检测 → 子载波间隔估计
[period_idx, delta_f_hat, detInfo] = detect_period_and_deltaf(R_use, df_psd, ...
    'FigureTitle', dataFile + " - 周期检测", 'PlotFigure', cfg.plotFigures, ...
    'MinPeakProminence', cfg.minPeakProminence, 'SmoothingSpan', cfg.smoothingSpan);

%% 6. 子载波数估计
if isnan(delta_f_hat)
    N_hat = NaN;
    subcarrierInfo = struct('continuous', NaN, 'adjusted', NaN, 'guard', NaN);
else
    [N_hat, subcarrierInfo] = estimate_subcarrier_number(B_hat, delta_f_hat);
end

isOFDM = detInfo.periodicityMetric >= periodicityThreshold && ~isnan(delta_f_hat);

%% 7. 汇总输出
fprintf('\n===== 估计结果 =====\n');
fprintf('数据文件: %s\n', dataFile);
fprintf('采样率 fs = %.3f MHz (请确认与实际一致)\n', fs/1e6);
fprintf('再处理方式: %s\n', reprocessDesc);
fprintf('估计带宽 B_hat = %.3f MHz (%.3f ~ %.3f MHz)\n', B_hat/1e6, f_L/1e6, f_H/1e6);
fprintf('周期索引 M_B = %.1f 样本 -> Δf_hat = %.3f kHz\n', period_idx, delta_f_hat/1e3);
fprintf('子载波数估计 N_hat = %.1f → 四舍五入 %d (_guard=%d)\n', ...
    subcarrierInfo.continuous, N_hat, subcarrierInfo.guard);
fprintf('周期性指标 = %.2f (阈值 %.2f) -> 判决: %s\n', detInfo.periodicityMetric, ...
    periodicityThreshold, ternary_str(isOFDM,'存在OFDM结构','未检测到OFDM结构'));

%% 8. 可视化汇总
if cfg.plotFigures
    figure('Name','检测流程总览');
    subplot(3,1,1);
    plot(f_grid/1e6, 10*log10(Pxx+eps)); grid on;
    xlabel('频率 (MHz)'); ylabel('PSD (dB/Hz)');
    title('Welch功率谱与带宽'); hold on;
    yl = ylim;
    plot([f_L f_L]/1e6, yl, 'r--');
    plot([f_H f_H]/1e6, yl, 'r--');
    legend('P_{xx}', 'f_L','f_H');

    subplot(3,1,2);
    h_basic = plot(axis_basic, R_basic, 'Color',[0.6 0.6 0.6]); hold on; grid on;
    h_selected = plot(R_use, 'b', 'LineWidth', 1.1);
    xlabel('索引 k'); ylabel('|G(k)|^2 (归一化)');
    title(['功率谱再处理对比 - ' reprocessDesc]);
    if numel(h_basic) > 1, h_basic = h_basic(1); end
    if numel(h_selected) > 1, h_selected = h_selected(1); end
    legend([h_basic, h_selected], {'basic','selected'});

    subplot(3,1,3);
    plot(detInfo.envelope, 'b'); hold on; grid on;
    if ~isempty(detInfo.peaks) && all(detInfo.peaks<=numel(detInfo.envelope))
        stem(detInfo.peaks, detInfo.envelope(detInfo.peaks), 'r');
    end
    xlabel('索引'); ylabel('平滑包络');
    title(sprintf('周期检测: M_B=%.1f, Δf=%.3f kHz, 指标=%.2f', ...
        period_idx, delta_f_hat/1e3, detInfo.periodicityMetric));
end

fprintf('\n完成。如需分析其他文件，请在 run.m 中修改 cfg。\n');

results = struct('dataFile', dataFile, 'fs', fs, 'B_hat', B_hat, 'f_L', f_L, ...
    'f_H', f_H, 'delta_f_hat', delta_f_hat, 'N_hat', N_hat, 'isOFDM', isOFDM, ...
    'detInfo', detInfo, 'subcarrierInfo', subcarrierInfo, 'cfg', cfg, ...
    'reprocessDesc', reprocessDesc, 'period_idx', period_idx);

end

%% 辅助: 默认参数
function cfg = apply_main_defaults(cfg)
    defaults = struct('dataFile','test3.dat', ...
        'fs',15e6, ...
        'maxSamples',2^21, ...
        'welchLen',4096, ...
        'welchWindow','hann', ...
        'welchOverlap',0.5, ...
        'useImprovedReprocess',true, ...
        'alphaIndex',1, ...
        'periodicityThreshold',2.0, ...
        'bandwidthThreshold_dB',6, ...
        'bandwidthGuardFraction',0.08, ...
        'bandwidthMinWidthHz',[], ...
        'bandwidthEnergyCoverage',0.995, ...
        'bandwidthDetectionMode','hybrid', ...
        'plotFigures',true, ...
        'zeroPadFactor',2, ...
    'suppressZeroLagPeak',true, ...
        'minPeakProminence',0.02, ...
        'smoothingSpan',51, ...
        'closeAllFigures',true, ...
        'clearCommandWindow',true);

    fn = fieldnames(defaults);
    for k = 1:numel(fn)
        name = fn{k};
        if ~isfield(cfg, name) || isempty(cfg.(name))
            cfg.(name) = defaults.(name);
        end
    end
end

%% 辅助: 打印友好字符串
function out = ternary_str(cond, yesText, noText)
    if cond
        out = yesText;
    else
        out = noText;
    end
end
