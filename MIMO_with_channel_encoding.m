%lmmse����㷨����ŵ��������ܣ����Ʒ�ʽΪQPSK���޿�ʱ����
clear all
Nt=4; %����������-----------2/4/8����������������������������
Nr=4; %����������----------------2/4/8��������������������������
msg=6;%ÿ����Ϣ���ݳ��ȣ��������������Ϊframe����
frame = 10; %ÿ֡�ĳ���
L = 10000; %�������֡��
EbN0 = 0:2:20;
M = 4; %QPSK����
a = randi([0,1],msg*L,Nt);%��Դ����
x = zeros(frame*L,Nt); %�ŵ����������ݣ��þ��󱣴�

for id=1:Nt
    for id2=1:L
        tempy=a((id2-1)*msg+1:id2*msg,id);
        tempy2=hamming_encode(tempy');
        x((id2-1)*frame+1:id2*frame,id)=tempy2';
    end
end
s = qpskmod_mine(x,M); %QPSK���ƺ���ź�
 %�����ŵ������໥���ŵ�����
 
%   testtest = randi([0,1],1,6);%���Ա��볤��
%   test_encode=hamming_encode(testtest);
%   test_decode=hamming_decode(test_encode);
%   [temp,ber0] = biterr(testtest,test_decode);
%   disp(ber0);
 
for index=1:length(EbN0)
    x_lmmse = [];
    for index1 = 1:L
        h = (randn(Nt,Nr)+j*randn(Nt,Nr)./sqrt(2)); %Rayleigh˥���ŵ���ÿһ֡����һ�β���
        sigma = sqrt(1/(10.^(EbN0(index)/10))); %ÿ���������ߵĸ�˹��������׼��
        n = sigma*(randn(frame,Nr)+j*randn(frame,Nr)); %ÿ���������ߵĸ�˹������
        w_lmmse = h'*inv(h*h'+1/(10.^(EbN0(index)/10))*diag(ones(1,Nt)));%w�����Ž⣨lmmse��,b=1 for qpsk
        y = s((index1-1)*frame+1:index1*frame,:)*h+n; %�ź�ͨ���ŵ�����ÿ�����ߵĽ��շֱ����
        y_lmmese =  y*w_lmmse; %�޸�������ʱ��lMMSE���
        temp_lmmse = qpskdemod_mine(y_lmmese,M); %�޸�������ʱ�Ľ����ÿ�������ϵĽ��
        x_lmmse = [x_lmmse;temp_lmmse]; %�޸�������ʱ�Ľ����� 
    end
    special_index1 = find(x_lmmse==2);x_lmmse(special_index1) = 0;%���ȷ�Χʱ������
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

figure%����ͼ��
semilogy(EbN0,ber,'-ro'); 
grid on;
legend('LMMSE');
title('4*4MIMO���ýṹ��Ϻ����ŵ�����LMMSE����㷨�Ĳ�������');%---------%22,44����һ��
xlabel('�����Eb/N0(dB)');
ylabel('������ʣ�BER��');



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