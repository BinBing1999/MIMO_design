%MMSE��lmmse����㷨���ܣ����Ʒ�ʽΪQPSK��v-blast��ϵ�����ļ����ڹ۲�����������Ӱ��
clear all
Nt = 2; %����������-----------2/4/8����������������������������
Nr = 2; %����������----------------2/4/8��������������������������
N = 10; %ÿ֡�ĳ���
L = 10000; %�������֡��
EbN0 = 0:2:20;
M = 4; %QPSK����
x = randi([0,1],N*L,Nt); %��Դ����
s = qpskmod_mine(x,M); %QPSK���ƺ���ź�
 %�����ŵ������໥���ŵ�����
for index=1:length(EbN0)
    x_mmse = [];
    x_lmmse = [];
    for index1 = 1:L
        h = (randn(Nt,Nr)+j*randn(Nt,Nr)./sqrt(2)); %Rayleigh˥���ŵ���ÿһ֡����һ�β���
        
        sigma1 = sqrt(1/(10.^(EbN0(index)/10))); %ÿ���������ߵĸ�˹��������׼��
        n = sigma1*(randn(N,Nr)+j*randn(N,Nr)); %ÿ���������ߵĸ�˹������
        w_mmse = h'*inv(h*h'+2*sigma1.^2*diag(ones(1,Nt))); %w�����Ž⣨mmse��
        w_lmmse = h'*inv(h*h'+1/(10.^(EbN0(index)/10))*diag(ones(1,Nt)));%w�����Ž⣨lmmse��,b=1 for qpsk
        y = s((index1-1)*N+1:index1*N,:)*h+n; %�ź�ͨ���ŵ�����ÿ�����ߵĽ��շֱ����
        
        
        y_mmese =  y*w_mmse; %�޸�������ʱ��MMSE���
        temp_mmse = qpskdemod_mine(y_mmese,M); %�����ÿ�������ϵĽ��
        x_mmse = [x_mmse;temp_mmse]; %������
        y_lmmese =  y*w_lmmse; %�޸�������ʱ��LMMSE���
        temp_lmmse = qpskdemod_mine(y_lmmese,M); %�����ÿ�������ϵĽ��
        x_lmmse = [x_lmmse;temp_lmmse]; %������
        
    end
    [temp,ber1(index)] = biterr(x,x_mmse,log2(M)); 
    [temp,ber2(index)] = biterr(x,x_lmmse,log2(M)); 
end

figure%����ͼ��
semilogy(EbN0,ber1,'-ko',EbN0,ber2,'-ro'); 
grid on;
legend('mmse','lmmse');
title('2*2MIMO���ýṹ��ͬ����㷨����');
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
