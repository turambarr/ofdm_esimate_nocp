function [bandwidth, f_low, f_high] = estimate_bandwidth(freq, Pxx, threshold_ratio)
% estimate_bandwidth  根据功率谱估计信号带宽及上下限
%
%   [bandwidth, f_low, f_high] = estimate_bandwidth(freq, Pxx, threshold_ratio)
%
% 输入参数:
%   freq           - 频率轴向量，与功率谱Pxx具有相同长度，单位Hz。
%   Pxx            - 平均功率谱密度向量。
%   threshold_ratio- 阈值比例，用于选取信号能量范围。例如0.1表示取最大值10%的阈值。
%
% 输出参数:
%   bandwidth      - 估计的信号带宽，单位Hz。
%   f_low          - 带宽下端频率，单位Hz。
%   f_high         - 带宽上端频率，单位Hz。
%
% 方法:
%   本函数通过对功率谱中能量强于某阈值的频率范围进行截取，估计信号的有效带宽。
%   具体而言，首先将功率谱归一化，根据最大值乘以threshold_ratio得到阈值，
%   然后找到功率谱大于阈值的最小索引和最大索引对应的频率即为信号频带范围。
%   最终带宽为频率差值(f_high - f_low)。

% 判断阈值比例是否合理
if nargin < 3 || isempty(threshold_ratio)
    threshold_ratio = 0.1; % 默认阈值为最大值的10%
end

% 归一化功率谱并计算阈值
P_norm = Pxx ./ max(Pxx);
threshold = threshold_ratio;

% 找到超过阈值的索引范围
indices = find(P_norm >= threshold);
if isempty(indices)
    % 如果没有找到，说明阈值过高，直接将带宽设为全频带
    f_low = freq(1);
    f_high = freq(end);
else
    f_low = freq(indices(1));
    f_high = freq(indices(end));
end

% 计算带宽
bandwidth = f_high - f_low;

end