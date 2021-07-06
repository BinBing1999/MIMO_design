%Alamouti 2��2�� ��ʱ�������ܣ����Ʒ�ʽΪQPSK
clear all
datasize = 100000; %����ķ�����
EbN0 = 0:2:20; %�����
M = 4; %QPSK����
x = randsrc(2,datasize/2,[0:3]); %����Դ����
x1 = qpskmod_mine(x,M);
h = randn(4,datasize/2) +j*randn(4,datasize/2); %Rayleigh˥���ŵ�
h = h./sqrt(2); %������һ��
ber=zeros(1,10);
for index=1:length(EbN0)
    sigma = sqrt(1/(2*10.^(EbN0(index)/10))); %Alamouti����ÿ�����ŵ���˹��������׼��
    n = sigma*(randn(4,datasize/2)+j*randn(4,datasize/2));
    n1(1,:) = (conj(h(1,:)).*n(1,:)+h(2,:).*conj(n(2,:))+conj(h(3,:)).*n(3,:)+h(4,:).*conj(n(4,:)))./(sum(abs(h).^2)); %2��2��Alamouti�����о�������������
    n1(2,:) = (conj(h(2,:)).*n(1,:)-h(1,:).*conj(n(2,:))+conj(h(4,:)).*n(3,:)-h(3,:).*conj(n(4,:)))./(sum(abs(h).^2));
    y1= x1 + n1;
    x_ans = qpskdemod_mine(y1,M);
    [temp,ber(index)] = biterr(x,x_ans,log2(M));
end
 
semilogy(EbN0,ber,'-bo')
grid on
legend('2��2��Alamouti')
xlabel('�����Eb/N0')
ylabel('������ʣ�BER��')
title('2��2��Alamouti���������˥�估��˹�������ŵ��µ�����')


%qpsk���ƺ���
function y = qpskmod_mine(x,M)
sita = 2*pi*x/M;
y = exp(1i*sita);
if isreal(y)
    y = complex(y,0);
end
end
%qpsk�������
function z = qpskdemod_mine(y,M)
y = y .* exp(1);
temp = M/(2*pi); 
z = round((angle(y) .* temp));
z(z < 0) = M + z(z < 0);
end

%�������룬dataÿ��һ�β���64λ��data
function hamming_code = hamming_encode(data)
% ��Ϣλ����
data_len = length(data);
% ��ලλ����flag_len
flag_len = 0;
while 2^flag_len - flag_len - 1 < data_len
    flag_len = flag_len + 1;
    flag_index(flag_len) = 2^(flag_len - 1);
end
% ���볤
code_len = data_len + flag_len;
% ��ʼ����Ԫ����
hamming_code = zeros(1, code_len);
% ����Դ��Ԫ
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

% �жϼලλ��ֵ
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
% ����Ϣλ���Ⱥͼලλ����
flag_len = 0;
while 2^flag_len < code_len
    flag_len = flag_len + 1;
end
data_len = code_len - flag_len;
% �ж���żУ��λ������ż�Բ�����
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
% ��λ������Ԫλ�ò�����
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