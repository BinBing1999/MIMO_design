%Alamouti 2发2收 空时编码性能，调制方式为QPSK
clear all
datasize = 100000; %仿真的符号数
EbN0 = 0:2:20; %信噪比
M = 4; %QPSK调制
x = randsrc(2,datasize/2,[0:3]); %数据源符号
x1 = qpskmod_mine(x,M);
h = randn(4,datasize/2) +j*randn(4,datasize/2); %Rayleigh衰落信道
h = h./sqrt(2); %能量归一化
ber=zeros(1,10);
for index=1:length(EbN0)
    sigma = sqrt(1/(2*10.^(EbN0(index)/10))); %Alamouti方案每个子信道高斯白噪声标准差
    n = sigma*(randn(4,datasize/2)+j*randn(4,datasize/2));
    n1(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+conj(h(3,:)).*n(3,:)+h(4,:).*conj(n(4,:)))./(sum(abs(h).^2)); %2发2收Alamouti方案判决变量的噪声项
    n1(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:))+conj(h(4,:)).*n(3,:)-h(3,:).*conj(n(4,:)))./(sum(abs(h).^2));
    y1= x1 + n1;
    x_ans = qpskdemod_mine(y1,M);
    [temp,ber(index)] = biterr(x,x_ans,log2(M));
end
 
semilogy(EbN0,ber,'-bo')
grid on
legend('2发2收Alamouti')
xlabel('信噪比Eb/N0')
ylabel('误比特率（BER）')
title('2发2收Alamouti设计在瑞利衰落及高斯白噪音信道下的性能')


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