%stbc of 4*4 mimo һ�������������������ʱ���룬��ʽ���ĵ�����alamouti���
clear all
M=4;
datasize = 999999; % ����ķ�����
EbNo = 0:2:20; % ����Ⱥ�����
x = randsrc(3,datasize/3,[0:3]); % 0-3�����
mm = qpskmod_mine(x,M);% QPSK����
ber=zeros(1,10);
h=randn(16,datasize/3)+j*randn(16,datasize/3); % ����˥���ŵ�˥��ϵ��sigma=1
h=h./sqrt(2); 
 
for index=1:length(EbNo)
    sigma=sqrt(1/(4/3*10.^(EbNo(index)/10))); % ÿ�����ŵ���˹��������׼�ȡ���������������
    n=sigma*(randn(16,datasize/3)+j*randn(16,datasize/3));

    %����Ϊ4*4mimo �ض���ʱ����
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
xlabel('����� EbNo(dB)')
ylabel('������� (BER)')
title('4*4����3ʱ������ʱ������������˥�估��˹�������ŵ��µ����� ')

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