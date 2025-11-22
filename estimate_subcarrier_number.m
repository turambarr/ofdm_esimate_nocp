function [N_hat, info] = estimate_subcarrier_number(B_hat, delta_f_hat, varargin)
% estimate_subcarrier_number  根据信号带宽与子载波间隔估计子载波数
%
%   [N_hat, info] = estimate_subcarrier_number(B_hat, delta_f_hat, Name, Value)
%
% 输入参数:
%   B_hat       - 估计的信号带宽 (Hz)。
%   delta_f_hat - 估计的子载波间隔 (Hz)。
%
% Name-Value参数:
%   'GuardSubcarriers' - 保护或空载波数量，默认2。
%
% 输出参数:
%   N_hat - 估计的子载波数，四舍五入取整。
%   info  - 中间信息结构体。
%
% 理论说明:
%   理想OFDM信号满足 B \approx N * Δf。盲估计场景下使用 round(B/Δf)
%   能在均方意义下最小化误差；引入 GuardSubcarriers 以扣除边缘保留的
%   空子载波，可减少对侧瓣/多径的敏感度。

validateattributes(B_hat, {'double','single'}, {'scalar','nonnegative'}, mfilename, 'B_hat');
validateattributes(delta_f_hat, {'double','single'}, {'scalar','positive'}, mfilename, 'delta_f_hat');

p = inputParser;
p.addParameter('GuardSubcarriers', 2, @(v) isnumeric(v) && isscalar(v) && v >= 0);
p.parse(varargin{:});
opts = p.Results;

continuous_est = B_hat / delta_f_hat;
adjusted_est = max(continuous_est - opts.GuardSubcarriers, 1);
N_hat = max(1, round(adjusted_est));

info = struct('continuous', continuous_est, ...
              'adjusted', adjusted_est, ...
              'guard', opts.GuardSubcarriers);

end