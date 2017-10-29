
% fp1 = fopen('D:\DWD\data\colDataout2.bin','rb');
% data3 = fread(fp1,8192*2,'float32');
% fp2 = fopen('D:\DWD\data\colDatain2.bin','rb');
% data4 = fread(fp2,8192*2,'float32');
% a = data4-data3;
% max(a)
% plot(a);
% figure;plot(data3);
j00=sqrt(-1);        % square root of -1
fp1 = fopen('D:\2016_8_DWD\20170923\data\g_tr1.bin','rb');
data3 = fread(fp1,16384,'float');
% fp2 = fopen('D:\2016_12_0507\DBF_Para_addr761_0509.bin','rb');
% data4 = fread(fp2,61*2,'float32');
a = data3(1:2:end) + j00*data3(2:2:end);
% b = data4(1:2:end) + j00*data4(2:2:end);
% b = fft(a);
% abs_b = abs(b);
% figure;plot(data3,data4,'b','r');
% plot(data4,'r');
% max(abs_b)
% figure;plot(abs(a));
% fp2 = fopen('D:\DWD\data\echobackup122612.bin','rb');
% data4 = fread(fp2,32768,'int8');
% a1 = data4(1:2:end) + j00*data4(2:2:end);
% figure;plot(abs(a1))
% max = data3 - data4
figure;plot(data3);





