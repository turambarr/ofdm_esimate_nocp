function [B_hat, f_L, f_H, info] = estimate_bandwidth_from_psd(Pxx, f_grid, varargin)
% estimate_bandwidth_from_psd  根据功率谱的能量分布自适应估计信号带宽
%
%   [B_hat, f_L, f_H, info] = estimate_bandwidth_from_psd(Pxx, f_grid, Name, Value)
%
% 输入参数:
%   Pxx    - 一维功率谱密度 (线性刻度)。
%   f_grid - 与Pxx对应的频率轴 (Hz)，长度需一致。
%
% Name-Value可选参数:
%   'GuardFraction'    - 估计噪声底时用的两侧保护带比例，默认0.1。
%   'Threshold_dB'     - 相对噪声底的门限 (dB)，默认6 dB。
%   'MinWidthHz'       - 最小带宽约束，默认 Δf。
%   'EnergyCoverage'   - 使用能量百分位(0~1)进一步收窄带宽，默认0.995。
%   'DetectionMode'    - 'threshold' / 'energy' / 'hybrid'，默认hybrid。
%   'PlotFigure'       - 是否绘图展示结果，默认true。
%   'FigureTitle'      - 绘图标题文本。
%
% 输出参数:
%   B_hat  - 估计的有效带宽 (Hz)。
%   f_L    - 低端频率 (Hz)。
%   f_H    - 高端频率 (Hz)。
%   info   - 结构体，包含噪声估计、门限、索引等中间变量。
%
% 背景:
%   论文提出利用积分功率谱寻找能量显著高于噪声的区间以确定带宽。本实现
%   先用谱两端区域估计噪声底Pn，再通过"Pn + α dB"的门限选出主瓣区间。
%   设α=6 dB时可在低SNR下保持鲁棒。若未检测到显著区间则退化为全频带。

validateattributes(Pxx, {'double','single'}, {'vector','real','nonnegative'}, mfilename, 'Pxx');
validateattributes(f_grid, {'double','single'}, {'vector'}, mfilename, 'f_grid');

if numel(Pxx) ~= numel(f_grid)
    error('estimate_bandwidth_from_psd:SizeMismatch', 'Pxx与f_grid长度不一致。');
end

Pxx = double(Pxx(:));
f_grid = double(f_grid(:));

df = mean(diff(f_grid));
if isnan(df) || df == 0
    df = 1;
end

p = inputParser;
p.addParameter('GuardFraction', 0.1, @(v) isnumeric(v) && v > 0 && v < 0.49);
p.addParameter('Threshold_dB', 6, @(v) isnumeric(v) && isscalar(v));
p.addParameter('MinWidthHz', abs(df), @(v) isnumeric(v) && isscalar(v) && v >= 0);
p.addParameter('PlotFigure', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('FigureTitle', '带宽估计', @(s) ischar(s) || isstring(s));
p.addParameter('EnergyCoverage', 0.995, @(v) isnumeric(v) && isscalar(v));
p.addParameter('DetectionMode', 'hybrid', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opts = p.Results;

N = numel(Pxx);
numGuard = max(1, round(N * opts.GuardFraction));
noiseIndices = [1:numGuard, N-numGuard+1:N];
Pn = median(Pxx(noiseIndices));
threshold = Pn * 10^(opts.Threshold_dB/10);

mask = Pxx >= threshold;
idx = find(mask);
if isempty(idx)
    f_L_thresh = f_grid(1);
    f_H_thresh = f_grid(end);
else
    f_L_thresh = f_grid(idx(1));
    f_H_thresh = f_grid(idx(end));
end

coverage = opts.EnergyCoverage;
if ~(coverage > 0 && coverage < 1)
    coverage = NaN;
end

energy = max(Pxx - Pn, 0);
if ~any(energy)
    energy = Pxx;
end
cumEnergy = cumsum(energy);
totalEnergy = cumEnergy(end);
energyAvailable = ~isnan(coverage) && totalEnergy > 0;
f_L_energy = NaN; f_H_energy = NaN;
if energyAvailable
    tailFrac = (1 - coverage) / 2;
    lowTarget = totalEnergy * tailFrac;
    highTarget = totalEnergy * (1 - tailFrac);
    idxLow = find(cumEnergy >= lowTarget, 1, 'first');
    idxHigh = find(cumEnergy >= highTarget, 1, 'first');
    if ~isempty(idxLow) && ~isempty(idxHigh)
        f_L_energy = f_grid(idxLow);
        f_H_energy = f_grid(idxHigh);
    else
        energyAvailable = false;
    end
end

mode = lower(string(opts.DetectionMode));
validModes = ["threshold","energy","hybrid"];
if ~any(mode == validModes)
    warning('estimate_bandwidth_from_psd:Mode', '未知DetectionMode=%s，已退回hybrid。', mode);
    mode = "hybrid";
end

switch mode
    case "threshold"
        f_L = f_L_thresh;
        f_H = f_H_thresh;
    case "energy"
        if energyAvailable
            f_L = f_L_energy;
            f_H = f_H_energy;
        else
            warning('estimate_bandwidth_from_psd:EnergyUnavailable', ...
                '能量百分位方法不可用，改用阈值法。');
            f_L = f_L_thresh;
            f_H = f_H_thresh;
        end
    otherwise % hybrid
        if energyAvailable
            f_cand_L = max(f_L_thresh, f_L_energy);
            f_cand_H = min(f_H_thresh, f_H_energy);
            if f_cand_H > f_cand_L
                f_L = f_cand_L;
                f_H = f_cand_H;
            else
                f_L = f_L_energy;
                f_H = f_H_energy;
            end
        else
            f_L = f_L_thresh;
            f_H = f_H_thresh;
        end
end

if (f_H - f_L) < opts.MinWidthHz
    halfWidth = opts.MinWidthHz/2;
    f_center = (f_H + f_L)/2;
    f_L = f_center - halfWidth;
    f_H = f_center + halfWidth;
end

B_hat = max(opts.MinWidthHz, f_H - f_L);

info = struct('noiseFloor', Pn, ...
             'threshold', threshold, ...
             'mask', mask, ...
             'indices', idx, ...
             'index_L_threshold', ternary(idx, @(v)v(1), NaN), ...
             'index_H_threshold', ternary(idx, @(v)v(end), NaN), ...
             'delta_f', df, ...
             'f_L_threshold', f_L_thresh, ...
             'f_H_threshold', f_H_thresh, ...
             'f_L_energy', f_L_energy, ...
             'f_H_energy', f_H_energy, ...
             'energyCoverage', coverage, ...
             'mode', mode);

if opts.PlotFigure
    figure('Name', char(opts.FigureTitle));
    plot(f_grid/1e6, 10*log10(Pxx + eps), 'LineWidth', 1.2); hold on; grid on;
    yl = ylim;
    plot([f_L f_L]/1e6, yl, 'r--', 'LineWidth', 1);
    plot([f_H f_H]/1e6, yl, 'r--', 'LineWidth', 1);
    plot(f_grid/1e6, 10*log10(info.noiseFloor)*ones(size(f_grid)), 'k-.');
    title(sprintf('%s | B_{hat}=%.2f MHz', opts.FigureTitle, B_hat/1e6));
    xlabel('频率 (MHz)'); ylabel('功率谱密度 (dB/Hz)');
    legend('P_{xx}', 'f_L', 'f_H', '噪声底');
info = add_energy_indices(info, cumEnergy, idxLow, idxHigh);
end

function out = ternary(cond, ifTrue, ifFalse)
    if isempty(cond)
        out = ifFalse;
    else
        out = ifTrue(cond);
    end
end

function info = add_energy_indices(info, cumEnergy, idxLow, idxHigh)
    if ~isempty(cumEnergy)
        info.index_L_energy = idxLow;
        info.index_H_energy = idxHigh;
    else
        info.index_L_energy = NaN;
        info.index_H_energy = NaN;
    end
end
end
