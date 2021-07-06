%lmmse检测算法结合信道编码性能，调制方式为QPSK，无空时编码
clear all
Nt=4; %发射天线数-----------2/4/8！！！！！！！！！！！！！！
Nr=4; %接收天线数----------------2/4/8！！！！！！！！！！！！！
msg=6;%每个信息内容长度，经过汉明编码变为frame长度
frame = 10; %每帧的长度
L = 10000; %仿真的总帧数
EbN0 = 0:2:20;
M = 4; %QPSK调制
a = randi([0,1],msg*L,Nt);%信源数据
x = zeros(frame*L,Nt); %信道编码后的数据，用矩阵保存

for id=1:Nt
    for id2=1:L
        tempy=a((id2-1)*msg+1:id2*msg,id);
        tempy2=hamming_encode(tempy');
        x((id2-1)*frame+1:id2*frame,id)=tempy2';
    end
end
s = qpskmod_mine(x,M); %QPSK调制后的信号
 %各个信道存在相互干扰的问题
 
%   testtest = randi([0,1],1,6);%测试编码长度
%   test_encode=hamming_encode(testtest);
%   test_decode=hamming_decode(test_encode);
%   [temp,ber0] = biterr(testtest,test_decode);
%   disp(ber0);
 
for index=1:length(EbN0)
    x_lmmse = [];
    for index1 = 1:L
        h = (randn(Nt,Nr)+j*randn(Nt,Nr)./sqrt(2)); %Rayleigh衰落信道，每一帧更新一次参数
        sigma = sqrt(1/(10.^(EbN0(index)/10))); %每根接收天线的高斯白噪声标准差
        n = sigma*(randn(frame,Nr)+j*randn(frame,Nr)); %每根接收天线的高斯白噪声
        w_lmmse = h'*inv(h*h'+1/(10.^(EbN0(index)/10))*diag(ones(1,Nt)));%w的最优解（lmmse）,b=1 for qpsk
        y = s((index1-1)*frame+1:index1*frame,:)*h+n; %信号通过信道，对每根天线的接收分别计算
        y_lmmese =  y*w_lmmse; %无干扰消除时的lMMSE检测
        temp_lmmse = qpskdemod_mine(y_lmmese,M); %无干扰消除时的解调，每根天线上的解调
        x_lmmse = [x_lmmse;temp_lmmse]; %无干扰消除时的解调结果 
    end
    special_index1 = find(x_lmmse==2);x_lmmse(special_index1) = 0;%过度范围时的修正
    special_index2 = find(x_lmmse==3);x_lmmse(special_index2) = 1;
    special_index3 = find(x_lmmse==4);x_lmmse(special_index3) = 0;
    special_index4 = find(x_lmmse==5);x_lmmse(special_index4) = 1;
    for id=1:Nt
        for id2=1:L
            tempy3=x_lmmse((id2-1)*frame+1:id2*frame,id);
            tempy4=hamming_decode(tempy3');
            a_lmmse((id2-1)*msg+1:id2*msg,id)=tempy4';
        end
    end
    [temp,ber(index)] = biterr(a,a_lmmse,log2(M)); 
end

figure%绘制图形
semilogy(EbN0,ber,'-ro'); 
grid on;
legend('LMMSE');
title('4*4MIMO复用结构结合汉明信道编码LMMSE检测算法的测试性能');%---------%22,44各测一次
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



%汉明编码，data每编一次产生64位新data
function hamming_code = hamming_encode(data)
% 信息位长度
data_len = length(data);
% 求监督位个数flag_len
flag_len = 0;
while 2^flag_len - flag_len - 1 < data_len
    flag_len = flag_len + 1;
    flag_index(flag_len) = 2^(flag_len - 1);
end
% 总码长
code_len = data_len + flag_len;
% 初始化码元数组
hamming_code = zeros(1, code_len);
% 填入源码元
index = 1;
i = 1;
while index <= code_len
    if ismember(index, flag_index)
        index = index + 1;
    else
        hamming_code(index) = data(i);
        index = index + 1;
        i = i + 1;
    end
end

% 判断监督位的值
for i = 1:flag_len
    flag_pos = 2^(i - 1);
    temp = 0;
    for j = flag_pos:2^i:code_len
        for index = j:min(j+2^(i-1)-1, code_len)
            temp = xor(temp, hamming_code(index));
        end
    end
    if temp == 1
        hamming_code(flag_pos) = 1;
    end
end
end

function output = hamming_decode(code)
code_len = length(code);
% 求信息位长度和监督位长度
flag_len = 0;
while 2^flag_len < code_len
    flag_len = flag_len + 1;
end
data_len = code_len - flag_len;
% 判断奇偶校验位组别的奇偶性并纠错
for i = 1:flag_len
    flag_pos = 2^(i - 1);
    temp = 0;
    for j = flag_pos:2^i:code_len
        for index = j:min(j+2^(i-1)-1, code_len)
            temp = xor(temp, code(index));
        end
    end
    if temp == 1
        wrong(i) = '1';
    else
        wrong(i) = '0';
    end
end
% 定位错误码元位置并纠错
wrong = reverse(wrong);
wrong_index = bin2dec(wrong);
if wrong_index >10
    wrong_index=0;
end
if wrong_index ~= 0
    code(wrong_index) = ~code(wrong_index);
end
for i = 1:flag_len
    flag_index(i) = 2^(i-1);
end
index = 0;
i = 1;
while index < code_len
    index = index + 1;
    if ismember(index, flag_index)
        continue
    else
        output(i) = code(index);
        i = i + 1;
    end
end
end