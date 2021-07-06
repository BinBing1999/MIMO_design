%stbc of 4*4 mimo 一种特殊的四天线正交空时编码，公式见文档，非alamouti设计
clear all
M=4;
datasize = 999999; % 仿真的符号数
EbNo = 0:2:20; % 信噪比横坐标
x = randsrc(3,datasize/3,[0:3]); % 0-3随机数
mm = qpskmod_mine(x,M);% QPSK调制
ber=zeros(1,10);
h=randn(16,datasize/3)+j*randn(16,datasize/3); % 瑞利衰落信道衰落系数sigma=1
h=h./sqrt(2); 
 
for index=1:length(EbNo)
    sigma=sqrt(1/(4/3*10.^(EbNo(index)/10))); % 每个子信道高斯白噪声标准差，取决于信噪比来计算
    n=sigma*(randn(16,datasize/3)+j*randn(16,datasize/3));

    %以下为4*4mimo 特定空时编码
    n1(1,:)=(conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+h(3,:).*conj(n(3,:))+conj(h(5,:)).*n(5,:)+h(6,:).*conj(n(6,:))+h(7,:).*conj(n(7,:))+conj(h(9,:)).*n(9,:)+h(10,:).*conj(n(10,:))+h(11,:).*conj(n(11,:))+conj(h(13,:)).*n(13,:)+h(14,:).*conj(n(14,:))+h(15,:).*conj(n(15,:)))./(sum(abs(h).^2)); 
    n1(2,:)=(-conj(h(1,:)).*n(2,:)+h(2,:).*conj(n(1,:))+h(4,:).*conj(n(3,:))-conj(h(5,:)).*n(6,:)+h(6,:).*conj(n(5,:))+h(8,:).*conj(n(7,:))-conj(h(9,:)).*n(10,:)+h(10,:).*conj(n(9,:))+h(12,:).*conj(n(11,:))-conj(h(13,:)).*n(14,:)+h(14,:).*conj(n(13,:))+h(16,:).*conj(n(15,:)))./(sum(abs(h).^2)); 
    n1(3,:)=(-conj(h(1,:)).*n(3,:)+h(3,:).*conj(n(1,:))-h(4,:).*conj(n(2,:))-conj(h(5,:)).*n(7,:)+h(7,:).*conj(n(5,:))-h(8,:).*conj(n(6,:))-conj(h(9,:)).*n(11,:)+h(11,:).*conj(n(9,:))-h(12,:).*conj(n(10,:))-conj(h(13,:)).*n(15,:)+h(15,:).*conj(n(13,:))-h(16,:).*conj(n(14,:)))./(sum(abs(h).^2)); 
    y = mm + n1;
    d = qpskdemod_mine(y,M);
    [temp,ber(index)]=biterr(x,d,log2(M));
end
 
semilogy(EbNo,ber,'-ko')
grid on
legend('QPSK')
xlabel('信噪比 EbNo(dB)')
ylabel('误比特率 (BER)')
title('4*4天线3时间间隔空时分组码在瑞利衰落及高斯白噪声信道下的性能 ')

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