function [t, ofdm_signal, params] = generate_ofdm_signal(num_subcarriers, num_symbols, cp_ratio, fs)
% generate_ofdm_signal  生成简单的OFDM基带信号用于算法验证
%
%   [t, ofdm_signal, params] = generate_ofdm_signal(num_subcarriers, num_symbols, cp_ratio, fs)
%
% 输入参数:
%   num_subcarriers - 子载波数目N。
%   num_symbols     - OFDM符号个数，用于生成较长的测试序列。
%   cp_ratio        - 循环前缀长度与有效符号长度的比例，例如0.25表示CP长度为符号长度的1/4。
%   fs              - 采样频率(Hz)。
%
% 输出参数:
%   t             - 时间向量(秒)，与生成的信号长度相同。
%   ofdm_signal    - 生成的复值基带OFDM信号。
%   params         - 结构体，包含子载波间隔delta_f、符号长度等信息。
%
% 方法:
%   本函数先生成随机QPSK调制的符号映射到N个子载波，通过IFFT得到时域OFDM符号。
%   然后在每个符号前添加循环前缀(CP)，将多个符号串联并输出连续时间信号。

% 计算每个OFDM符号的样本数
N = num_subcarriers;

% 有效符号长度为N个采样点，对应采样频率fs下的持续时间
Ts_symbol = N / fs;

% 循环前缀长度
cp_len = round(cp_ratio * N);

% 单个OFDM符号包含的总样本数
N_total = N + cp_len;

% 预分配输出向量
ofdm_signal = [];

% 子载波间隔
delta_f = fs / N;

for k = 1:num_symbols
    % 生成随机QPSK符号
    data_bits = randi([0 3], N, 1); % 四相位取值0~3
    % QPSK调制映射：0->1+1j, 1->1-1j, 2->-1+1j, 3->-1-1j
    qpsk_map = [1+1j, 1-1j, -1+1j, -1-1j] / sqrt(2);
    symbols = qpsk_map(data_bits + 1).';
    
    % 对每个OFDM符号进行IFFT
    time_symbol = ifft(symbols, N);
    
    % 添加循环前缀
    cp = time_symbol(end - cp_len + 1:end);
    ofdm_with_cp = [cp; time_symbol];
    
    % 累加到输出序列
    ofdm_signal = [ofdm_signal; ofdm_with_cp];
end

% 构造时间向量
num_samples = length(ofdm_signal);
t = (0:num_samples-1).' / fs;

% 填充参数结构体
params.delta_f = delta_f;
params.N = N;
params.cp_len = cp_len;
params.N_total = N_total;
params.Ts_symbol = Ts_symbol;

end