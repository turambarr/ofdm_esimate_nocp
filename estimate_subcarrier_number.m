function [N_est, delta_f] = estimate_subcarrier_number(bandwidth, fs, M)
% estimate_subcarrier_number  根据带宽和周期估计OFDM子载波数
%
%   [N_est, delta_f] = estimate_subcarrier_number(bandwidth, fs, M)
%
% 输入参数:
%   bandwidth - 估计的信号带宽 (Hz)。
%   fs        - 采样频率 (Hz)。
%   M         - “谱再处理”序列的周期，以样本点数表示。
%
% 输出参数:
%   N_est     - 估计的子载波数，四舍五入取整。
%   delta_f   - 估计的子载波间隔 (Hz)。
%
% 根据文献公式，谱再处理序列的主周期M对应的频率间隔为采样频率除以M，
% 即子载波间隔 delta_f = fs / M。信号带宽约等于子载波数乘以子载波间隔，
% 因此子载波数估计为带宽除以间隔，并取最近整数。

% 计算子载波间隔
delta_f = fs / M;

% 估计子载波数并取整
N_est = round(bandwidth / delta_f);

end