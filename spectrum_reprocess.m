function [reprocessed_spectrum] = spectrum_reprocess(Pxx)
% spectrum_reprocess  功率谱二次处理（功率谱的功率谱）
%
%   reprocessed_spectrum = spectrum_reprocess(Pxx)
%
% 输入参数:
%   Pxx  - 一维功率谱密度向量（单边或双边均可），长度为N。
%
% 输出参数:
%   reprocessed_spectrum - 经过二次处理后的功率谱，其实质为对输入功率谱
%                          进行FFT变换并取幅度平方，体现功率谱的周期性特征。
%
% 原理:
%   OFDM信号的功率谱在频域呈现准周期性，计算功率谱的FFT可以突出这种
%   周期结构。通过对功率谱进行FFT并求取幅度平方，即得到所谓的“谱再处理”。

% 对功率谱进行FFT运算
fft_P = fft(Pxx);

% 取幅度平方得到再处理后的谱
reprocessed_spectrum = abs(fft_P).^2;

end