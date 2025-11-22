function [R_basic, axis_idx] = spectrum_reprocessing_basic(Pxx, varargin)
% spectrum_reprocessing_basic  对功率谱执行FFT并取模平方以突出隐周期
%
%   [R_basic, axis_idx] = spectrum_reprocessing_basic(Pxx, Name, Value)
%
% 输入:
%   Pxx - 功率谱密度列向量。
%
% Name-Value参数:
%   'Normalize'    - 是否按最大值归一，默认true。
%   'ZeroPadFactor'- 再处理FFT的零填充倍数，默认2。
%   'PlotFigure'   - 是否绘图，默认true。
%   'FigureTitle'  - 绘图标题。
%
% 输出:
%   R_basic  - “功率谱的功率谱”，对应|FFT{Pxx}|^2。
%   axis_idx - 与R_basic对应的索引（延迟/周期样本）。
%
% 理论说明:
%   等间隔子载波使Pxx呈现准周期结构，相当于频域的梳状滤波器。对Pxx
%   再做一次FFT等价于计算其自相关，隐含的梳状周期在|FFT|^2后转为一列
%   均匀间隔的峰值。该处理对噪声ほ不放大，从而提高检测SNR。

validateattributes(Pxx, {'double','single'}, {'vector','real','nonnegative'}, mfilename, 'Pxx');
Pxx = double(Pxx(:));

p = inputParser;
p.addParameter('Normalize', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('ZeroPadFactor', 2, @(v) isnumeric(v) && v >= 1);
p.addParameter('PlotFigure', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('FigureTitle', '基本功率谱再处理', @(s) ischar(s) || isstring(s));
p.addParameter('RemoveMean', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('SuppressZeroLag', true, @(v) islogical(v) || isnumeric(v));
p.parse(varargin{:});
opts = p.Results;

N = numel(Pxx);
M = 2^nextpow2(N * opts.ZeroPadFactor);
if opts.RemoveMean
    proc = Pxx - mean(Pxx);
else
    proc = Pxx;
end
G = fft(proc, M);
R_basic = abs(G).^2;

if opts.SuppressZeroLag && ~isempty(R_basic)
    R_basic(1) = 0;          % k=0 (零延迟)
    R_basic(end) = 0;        % 等效的周期镜像
end

if opts.Normalize
    denom = max(R_basic + eps);
    if denom > 0
        R_basic = R_basic ./ denom;
    end
end

axis_idx = (0:M-1).';

if opts.PlotFigure
    figure('Name', char(opts.FigureTitle));
    R_plot = 10*log10(R_basic + eps);
    plot(axis_idx, R_plot, 'LineWidth', 1.2); grid on;
    xlabel('索引 k (对应虚拟延迟)'); ylabel('|G(k)|^2 (dB)');
    title(opts.FigureTitle + " - OFDM隐周期可视化");
end
end
