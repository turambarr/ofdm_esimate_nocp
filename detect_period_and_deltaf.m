function [period_idx, delta_f_hat, info] = detect_period_and_deltaf(R, freqResolution, varargin)
% detect_period_and_deltaf  在谱再处理结果中检测周期并估计子载波间隔
%
%   [period_idx, delta_f_hat, info] = detect_period_and_deltaf(R, freqResolution, Name, Value)
%
% 输入:
%   R              - 再处理后的谱 (basic或improved)，实非负向量。
%   freqResolution - 原始功率谱的频率分辨率 Δf_psd = fs / L。
%
% Name-Value参数:
%   'SmoothingSpan'     - 用于Envelope的移动平均窗口，默认31点。
%   'MinPeakProminence' - 峰值最小突出度，默认0.1。
%   'PlotFigure'        - 是否绘图展示峰检测，默认true。
%   'FigureTitle'       - 图标题。
%
% 输出:
%   period_idx   - 估计的周期 (单位:索引点)。
%   delta_f_hat  - 子载波间隔估计 (Hz)，满足 Δf_hat = Δf_psd / period_idx。
%   info         - 结构体，含峰位置、平滑包络、周期性指标等。
%
% 理论说明:
%   R(k)≈|\mathcal{F}{Pxx}|^2 时，其主峰间距与子载波间隔成倒数关系。
%   若功率谱分辨率为Δf_psd，则再处理域中的周期索引M_B与子载波间隔
%   满足 Δf_OFDM = Δf_psd / M_B。该函数通过平滑+峰值检测+自相关备选方案
%   获得M_B，并输出周期可信度供盲检测判决。

validateattributes(R, {'double','single'}, {'vector','real','nonnegative'}, mfilename, 'R');
validateattributes(freqResolution, {'double','single'}, {'scalar','positive'}, mfilename, 'freqResolution');

R = double(R(:));
R = R ./ max(R + eps);

p = inputParser;
p.addParameter('SmoothingSpan', 31, @(v) isnumeric(v) && isscalar(v) && v >= 5);
p.addParameter('MinPeakProminence', 0.1, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('PlotFigure', true, @(v) islogical(v) || isnumeric(v));
p.addParameter('FigureTitle', '周期检测', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opts = p.Results;

span = round(opts.SmoothingSpan);
span = span + mod(span+1,2); % 使用奇数窗口
R_env = movmean(R, span);
R_centered = R_env - mean(R_env);
R_positive = max(R_centered, 0);

[minProminence, minDist] = deal(opts.MinPeakProminence, max(3, round(span/2)));
[pkVals, pkLocs] = findpeaks(R_positive, 'MinPeakProminence', minProminence, ...
    'MinPeakDistance', minDist);

period_idx = NaN;
method = 'peaks';

if numel(pkLocs) >= 2
    intervals = diff(pkLocs);
    period_idx = round(mean(intervals));
else
    % 退化时采用频域法：对去直流的包络做FFT，寻找最大非零频率
    method = 'fft';
    spec = abs(fft(R_centered));
    spec(1) = 0;
    halfLen = floor(numel(spec)/2);
    [specMax, relIdx] = max(spec(2:halfLen));
    if ~isempty(relIdx) && specMax > 0
        bin = relIdx + 1;           % 对应原始FFT的频率索引
        harmonic = bin - 1;         % 相对零频的基波序号
        if harmonic >= 1
            period_idx = round(numel(R_centered) / harmonic);
            pkVals = spec(bin);
            pkLocs = bin;
        end
    end
end

if isnan(period_idx) || period_idx <= 0
    delta_f_hat = NaN;
else
    delta_f_hat = freqResolution / period_idx;
end

% 周期性指标: 主峰/均值
if isempty(pkVals)
    periodicity_metric = 0;
else
    periodicity_metric = max(pkVals) / max(std(R_centered), eps);
end

info = struct('peaks', pkLocs, ...
              'peakValues', pkVals, ...
              'intervals', exist_or_empty('intervals'), ...
              'envelope', R_env, ...
              'centered', R_centered, ...
              'method', method, ...
              'periodicityMetric', periodicity_metric);

if opts.PlotFigure
    figure('Name', char(opts.FigureTitle));
    subplot(2,1,1);
    hR = plot(R, 'Color', [0.6 0.6 0.6]); hold on;
    hEnv = plot(R_env, 'b', 'LineWidth', 1.2); grid on;
    if ~isempty(pkLocs)
        hPk = stem(pkLocs, R_env(pkLocs), 'r', 'filled');
        legend([hR hEnv hPk], {'R','平滑包络','峰'});
    else
        legend([hR hEnv], {'R','平滑包络'});
    end
    title(sprintf('%s | 周期候选 = %s', opts.FigureTitle, num2str(period_idx)));
    ylabel('归一化幅值'); xlabel('索引');

    subplot(2,1,2);
    acf_disp = ifft(abs(fft(R_env)).^2);
    acf_disp = real(acf_disp(:));
    acf_disp = acf_disp ./ max(acf_disp + eps);
    plot(acf_disp, 'LineWidth', 1.2); grid on;
    xlabel('索引'); ylabel('自相关 (归一化)');
    title(sprintf('辅助自相关 | 指标=%.2f', periodicity_metric));
end
end

function out = exist_or_empty(name)
% exist_or_empty  帮助函数，在变量尚未声明时返回[]
    if evalin('caller', sprintf('exist(''%s'', ''var'')', name))
        out = evalin('caller', name);
    else
        out = [];
    end
end
