function [freq, Pxx] = compute_power_spectrum(signal, fs, L)
% compute_power_spectrum  使用Welch法估计信号的平均功率谱
%
%   [freq, Pxx] = compute_power_spectrum(signal, fs, L)
%
% 输入参数:
%   signal - 输入的复数或实数信号向量（行或列向量均可），单位时间域采样值。
%   fs     - 采样频率（Hz）。
%   L      - Welch算法中每段的长度，用于控制频率分辨率，较大的L可以提高频率分辨率。
%
% 输出参数:
%   freq   - 与功率谱对应的频率轴（Hz），长度等于功率谱长度的一半加一，因为pwelch默认输出单边谱。
%   Pxx    - 由Welch法估计得到的平均功率谱密度 (功率/Hz)。
%
% 该函数利用MATLAB内置的pwelch函数进行功率谱密度估计，窗口函数使用汉宁窗，
% 为了提高频率分辨率通常选择适当的段长和重叠比例。pwelch返回的是单边谱，
% 因此频率范围为[0, fs/2]。

% 确保输入信号为列向量
signal = signal(:);

% 设置Welch方法参数: 使用汉宁窗、50%重叠
window = hanning(L);
nOverlap = floor(L/2);
nFFT = L;  % FFT点数与段长一致

% 调用pwelch估计功率谱密度
[Pxx, f] = pwelch(signal, window, nOverlap, nFFT, fs);

% 将频率轴输出
freq = f;

end