%MMSE与lmmse检测算法性能，调制方式为QPSK，v-blast体系，该文件用于观测天线数量的影响
clear all
Nt = 2; %发射天线数-----------2/4/8！！！！！！！！！！！！！！
Nr = 2; %接收天线数----------------2/4/8！！！！！！！！！！！！！
N = 10; %每帧的长度
L = 10000; %仿真的总帧数
EbN0 = 0:2:20;
M = 4; %QPSK调制
x = randi([0,1],N*L,Nt); %信源数据
s = qpskmod_mine(x,M); %QPSK调制后的信号
 %各个信道存在相互干扰的问题
for index=1:length(EbN0)
    x_mmse = [];
    x_lmmse = [];
    for index1 = 1:L
        h = (randn(Nt,Nr)+j*randn(Nt,Nr)./sqrt(2)); %Rayleigh衰落信道，每一帧更新一次参数
        
        sigma1 = sqrt(1/(10.^(EbN0(index)/10))); %每根接收天线的高斯白噪声标准差
        n = sigma1*(randn(N,Nr)+j*randn(N,Nr)); %每根接收天线的高斯白噪声
        w_mmse = h'*inv(h*h'+2*sigma1.^2*diag(ones(1,Nt))); %w的最优解（mmse）
        w_lmmse = h'*inv(h*h'+1/(10.^(EbN0(index)/10))*diag(ones(1,Nt)));%w的最优解（lmmse）,b=1 for qpsk
        y = s((index1-1)*N+1:index1*N,:)*h+n; %信号通过信道，对每根天线的接收分别计算
        
        
        y_mmese =  y*w_mmse; %无干扰消除时的MMSE检测
        temp_mmse = qpskdemod_mine(y_mmese,M); %解调，每根天线上的解调
        x_mmse = [x_mmse;temp_mmse]; %解调结果
        y_lmmese =  y*w_lmmse; %无干扰消除时的LMMSE检测
        temp_lmmse = qpskdemod_mine(y_lmmese,M); %解调，每根天线上的解调
        x_lmmse = [x_lmmse;temp_lmmse]; %解调结果
        
    end
    [temp,ber1(index)] = biterr(x,x_mmse,log2(M)); 
    [temp,ber2(index)] = biterr(x,x_lmmse,log2(M)); 
end

figure%绘制图形
semilogy(EbN0,ber1,'-ko',EbN0,ber2,'-ro'); 
grid on;
legend('mmse','lmmse');
title('2*2MIMO复用结构不同检测算法性能');
xlabel('信噪比Eb/N0(dB)');
ylabel('误比特率（BER）');



%qpsk调制函数
function y = qpskmod_mine(x,M)
sita = 2*pi*x/M;
y = exp(1i*sita);
if isreal(y)
    y = complex(y,0);
end
end
%qpsk解调函数
function z = qpskdemod_mine(y,M)
y = y .* exp(1);
temp = M/(2*pi); 
z = round((angle(y) .* temp));
z(z < 0) = M + z(z < 0);
end
