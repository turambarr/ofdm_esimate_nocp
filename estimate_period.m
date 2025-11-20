function [M, peak_locs] = estimate_period(reprocessed_spectrum)
% estimate_period  估计“谱再处理”序列的主周期
%
%   [M, peak_locs] = estimate_period(reprocessed_spectrum)
%
% 输入参数:
%   reprocessed_spectrum - “谱再处理”后的序列，是经过FFT和幅度平方得到的
%                          能量谱，一般为实数非负向量。
%
% 输出参数:
%   M         - 估计的周期（单位：点数），即相邻主峰之间的平均间隔。
%   peak_locs - 找到的峰值位置索引向量，用于辅助调试和绘图。
%
% 方法:
%   利用信号处理工具箱中的findpeaks函数检测能量谱中的峰值，限定峰值之间的
%   最小距离，得到明显的周期峰。计算相邻峰之间的距离并取平均值作为周期M。

% 使用findpeaks找到谱中的峰值
% 为避免检测过多杂峰，可以设定最低峰值高度为全局最大值的一定比例
peakThreshold = max(reprocessed_spectrum) * 0.5; % 50%阈值

% findpeaks需要Signal Processing Toolbox，在MATLAB标准发行版中一般包含
[pks, locs] = findpeaks(reprocessed_spectrum, 'MinPeakHeight', peakThreshold);

% 如果检测到的峰不足两个，则无法估计周期
if numel(locs) < 2
    M = NaN;
    peak_locs = locs;
    return;
end

% 计算相邻峰的间隔并求平均值
intervals = diff(locs);
M = round(mean(intervals));

peak_locs = locs;

end