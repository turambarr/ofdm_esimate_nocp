function [Pxx, f_grid, info] = estimate_psd_welch(x, fs, segmentLength, varargin)
% estimate_psd_welch  使用Welch方法估计功率谱并返回完整频率栅格
%
%   [Pxx, f_grid, info] = estimate_psd_welch(x, fs, segmentLength, Name, Value)
%
% 输入参数:
%   x             - 复基带或实信号列向量。
%   fs            - 采样率 (Hz)。
%   segmentLength - Welch算法中每段的采样点数 (即FFT长度 L)，
%                   决定频率分辨率 Δf = fs / L。
%
% Name-Value可选参数:
%   'Window'   - 窗函数类型，默认为 'hann'。
%   'Overlap'  - 段与段之间的重叠比例 (0~0.95)，默认 0.5。
%   'Nfft'     - FFT长度，默认等于 segmentLength。
%   'CenterDC' - 是否通过fftshift让0 Hz位于中心，默认 true。
%
% 输出参数:
%   Pxx    - Welch法得到的双边功率谱密度 (线性刻度)，长度等于Nfft。
%   f_grid - 与Pxx对应的频率轴 (Hz)。
%   info   - 结构体，包含:
%              .delta_f  频率分辨率 (Hz)
%              .window   实际使用的窗函数样本
%              .noverlap 实际的重叠点数
%              .nfft     使用的FFT长度
%
% 说明:
%   采用MATLAB内置pwelch，并使用'twosided'选项得到完整的双边谱，
%   再进行fftshift以便频率轴按[-fs/2, fs/2)排序。该实现作为后续
%   "功率谱再处理"模块的基础，需保持定长频率分辨率。
%
% 作者注:
%   OFDM盲检测依赖功率谱的隐周期性，本函数确保Pxx满足等距采样，
%   以便后续对Pxx再做FFT时仍然对应物理频率刻度。

p = inputParser;
p.addParameter('Window', 'hann', @(s) (ischar(s) || isstring(s)));
p.addParameter('Overlap', 0.5, @(v) isnumeric(v) && isscalar(v) && v >= 0 && v < 1);
p.addParameter('Nfft', segmentLength, @(v) isnumeric(v) && isscalar(v) && v > 0);
p.addParameter('CenterDC', true, @(v) islogical(v) || isnumeric(v));
p.parse(varargin{:});
opts = p.Results;


validateattributes(x, {'double','single'}, {'vector'}, mfilename, 'x');
validateattributes(fs, {'double','single'}, {'scalar','positive'}, mfilename, 'fs');
validateattributes(segmentLength, {'double','single'}, {'scalar','integer','>=',4}, mfilename, 'segmentLength');
x = double(x(:)); % 统一为列向量并强制双精度

sigLen = numel(x);
if segmentLength > sigLen
    segmentLength = max(32, 2^floor(log2(sigLen)));
    warning('estimate_psd_welch:SegmentLength', ...
        '段长大于信号长度，已调整为 %d 点。', segmentLength);
end

L = round(segmentLength);
nfft = max(round(opts.Nfft), L);
noverlap = min(round(L * opts.Overlap), L-1);
window = get_window(opts.Window, L);

% 使用Welch法获得双边功率谱
[P_raw, f_raw] = pwelch(x, window, noverlap, nfft, fs, 'twosided');

if opts.CenterDC
    Pxx = fftshift(P_raw);
    f_grid = (-nfft/2 : nfft/2 - 1).' * (fs / nfft);
else
    Pxx = P_raw;
    f_grid = f_raw;
end

info = struct('delta_f', fs / nfft, ...
              'window', window, ...
              'noverlap', noverlap, ...
              'nfft', nfft);
end

function w = get_window(name, L)
% get_window  根据名称生成窗函数，默认返回列向量。
    name = lower(string(name));
    switch name
        case {"hann","hanning"}
            w = hann(L, 'periodic');
        case "hamming"
            w = hamming(L, 'periodic');
        case "rect"
            w = rectwin(L);
        otherwise
            warning('estimate_psd_welch:get_window', ...
                '未知窗类型 %s，使用Hann窗替代。', name);
            w = hann(L, 'periodic');
    end
    w = w(:);
end
