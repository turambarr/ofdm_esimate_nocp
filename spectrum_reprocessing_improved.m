function R_improved = spectrum_reprocessing_improved(Pxx, varargin)
% spectrum_reprocessing_improved  提升多径鲁棒性的功率谱再处理变体
%
%   R_improved = spectrum_reprocessing_improved(Pxx, Name, Value)
%
% 该函数在基本谱再处理的基础上实现论文提出的两类改进:
%   A) 对功率谱取α次根(α=1/2或1/4)，抑制乘性通道波动再进行FFT。
%   B) 先对时域信号求解析信号的瞬时幅度，再计算其功率谱并执行改进A。
%
% Name-Value参数:
%   'Alphas'          - α列表 (默认[0.5 0.25])。
%   'ZeroPadFactor'   - FFT零填充倍数，默认2。
%   'Normalize'       - 是否归一化输出，默认true。
%   'PlotFigure'      - 是否绘制比较图，默认true。
%   'FigureTitle'     - 图标题。
%   'TimeSignal'      - 原始时域信号x[n]，启用方式B所必需。
%   'Fs'              - 采样率Hz，用于方式B的Welch参数。
%   'SegmentLength'   - 方式B使用的Welch段长，默认1024。
%   'Overlap'         - 方式B Welch重叠比例，默认0.5。
%
% 输出:
%   R_improved - 结构体，包含:
%                 .root{k}      : α_k对应的谱再处理
%                 .alphas       : 使用的α列表
%                 .instAmp      : 基于瞬时幅度的谱再处理(若可用)
%                 .axis         : 索引轴
%
% 理论背景:
%   多径信道使功率谱出现缓慢起伏，相当于Pxx乘以|H(f)|^2。对Pxx取α次根时，
%   利用(1+x)^α的泰勒展开可知乘性扰动x被压缩，从而减小对隐周期的掩盖。
%   解析信号瞬时幅度仅保留幅度包络，过滤频率调制与载频项，使重构的Pxx
%   主要反映子载波功率分布，进一步提升周期峰清晰度。

validateattributes(Pxx, {'double','single'}, {'vector','real','nonnegative'}, mfilename, 'Pxx');
Pxx = double(Pxx(:));

p = inputParser;
p.addParameter('Alphas', [0.5 0.25], @(v) isnumeric(v) && all(v>0) && all(v<=1));
p.addParameter('ZeroPadFactor', 2, @(v) isnumeric(v) && v>=1);
p.addParameter('Normalize', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('PlotFigure', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('FigureTitle', '改进功率谱再处理', @(s) ischar(s) || isstring(s));
p.addParameter('TimeSignal', [], @(v) isnumeric(v));
p.addParameter('Fs', [], @(v) isnumeric(v) && isscalar(v) && v>0);
p.addParameter('SegmentLength', 1024, @(v) isnumeric(v) && isscalar(v) && v>=32);
p.addParameter('Overlap', 0.5, @(v) isnumeric(v) && isscalar(v) && v>=0 && v<1);
p.addParameter('RemoveMean', true, @(v) islogical(v) || isnumeric(v));
p.parse(varargin{:});
opts = p.Results;

alphas = unique(opts.Alphas(:).');
N = numel(Pxx);
M = 2^nextpow2(N * opts.ZeroPadFactor);
axis_idx = (0:M-1).';

rootSpectra = cell(size(alphas));
for k = 1:numel(alphas)
    alpha = alphas(k);
    P_alpha = (Pxx + eps).^alpha; % 非线性压缩
    if opts.RemoveMean
        P_alpha = P_alpha - mean(P_alpha);
    end
    G_alpha = fft(P_alpha, M);
    R_alpha = abs(G_alpha).^2;
    if opts.Normalize
        R_alpha = R_alpha ./ max(R_alpha + eps);
    end
    rootSpectra{k} = R_alpha;
end

instAmpSpectrum = [];
if ~isempty(opts.TimeSignal) && ~isempty(opts.Fs)
    ts = opts.TimeSignal(:);
    if isreal(ts)
        xa = hilbert(ts);
    else
        xa = ts; % 已是解析信号
    end
    instAmp = abs(xa);
    [P_amp, ~, meta] = estimate_psd_welch(instAmp, opts.Fs, opts.SegmentLength, ...
        'Overlap', opts.Overlap, 'CenterDC', true);
    alpha_amp = min(alphas); % 使用最小α获得最强抑制
    P_amp_root = (P_amp + eps).^alpha_amp;
    if opts.RemoveMean
        P_amp_root = P_amp_root - mean(P_amp_root);
    end
    G_amp = fft(P_amp_root, M);
    instAmpSpectrum = abs(G_amp).^2;
    if opts.Normalize
        instAmpSpectrum = instAmpSpectrum ./ max(instAmpSpectrum + eps);
    end
else
    meta = struct('delta_f', NaN);
end

R_improved = struct('root', {rootSpectra}, ...
                    'alphas', alphas, ...
                    'instAmp', instAmpSpectrum, ...
                    'axis', axis_idx, ...
                    'delta_f_psd', meta.delta_f);

if opts.PlotFigure
    figure('Name', char(opts.FigureTitle));
    tiledlayout(numel(alphas) + (~isempty(instAmpSpectrum)), 1, 'TileSpacing','compact');
    for k = 1:numel(alphas)
        nexttile;
        plot(axis_idx, rootSpectra{k}, 'LineWidth', 1.2); grid on;
        title(sprintf('取根演化方案 (\alpha = %.2f)', alphas(k)));
        xlabel('索引'); ylabel('归一化功率');
    end
    if ~isempty(instAmpSpectrum)
        nexttile;
        plot(axis_idx, instAmpSpectrum, 'LineWidth', 1.2); grid on;
        title('解析信号瞬时幅度谱再处理');
        xlabel('索引'); ylabel('归一化功率');
    end
    sgtitle(opts.FigureTitle);
end
end
