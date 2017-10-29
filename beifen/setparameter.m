%��������
fd = fopen('E:\L_20170609\LSAR_Stripmap_Br200MHz_Tr20us_Fsr250MHz_PRF1000Hz_32K_20170609_1205_4_0.dat','rb');
% fd = fopen('I:\DWD\L��������\LSAR_Stripmap_Br200MHz_Tr20us_Fsr250MHz_PRF1000Hz_32K_20170609_1337_29_0.dat','rb');
% fseek(fd,(32768*4+128*8)*32768*2,'bof');

%% ������ͷ
Ass = fread(fd,128,'uint8');
index = Ass(28);
Tmin = ( Ass(47) + Ass(48)*256 )*1e-6;      %������ʱ
Tp = Ass(27,1)*1e-6;                        %������20us
fc = 1300*1e6;             %L������Ƶ
% fc = 390*1e6;               %P������Ƶ
Br = (Ass(34)*256+Ass(33))*1e6;              %����200M
Fs = (Ass(36)*256+Ass(35))*1e6;             %������250M
prf = Ass(32)*256+Ass(31);                  %PRF1000
Nr = (Ass(30)*256+Ass(29))*1024;            %�����������
Nr = Nr/2;
% Nr = 8192;
H = (Ass(66)*256^3+Ass(65)*256^2+Ass(64)*256+Ass(63))*0.1;
% angle_elevation = Ass(74)*256+Ass(73);      %�������ĸ�����
% if angle_elevation>0
%     angle_elevation = -1*angle_elevation*pi/180;
% end
% henggun = (Ass(50)*256+Ass(49))*0.001*pi/180;     %�ɻ������
% fy = (Ass(52)*256+Ass(51))*0.001*pi/180;          %�ɻ�������
% hangxiang = (Ass(54)*256+Ass(53))*0.01*pi/180;   %�ɻ������?
c = 3e8;                %����
j = sqrt(-1);           %j^2=-1

Rmin = Tmin*c/2;
% Rmin = 3600;
va = 400/3.6;           %�ٶ�  
D = 1.36;           %L�������߿ھ�
% D = 1.44*0.4;           %P�������߿ھ�

lamda = c/fc;           %����
bin_r = c/2/Fs;         %���������������
gama = Br/Tp;            %�������Ƶ��?
delta_theta = 0.89*lamda/D;  %��λ�������
delta1 = rad2deg(delta_theta);

R = Rmin+((1:Nr)-1)*bin_r;      %�����������þ���


% angle_azimuth = pi/2;        %�������ķ�λ��
% [ReYzPit,ReYzAzim] = ReAirSpace(angle_azimuth, angle_elevation, henggun, fy, hangxiang);            %ת���������ϵ��
% thetas_c = acos(cos(ReYzAzim)*cos(ReYzPit));            %ǰб��
% fdc = 2*va*cos(thetas_c)/lamda;                         %����������Ƶ��
% fdr = -2*va*va*sin(thetas_c)^2/lamda/Ra;                %�����յ�Ƶ��
% Bd = abs(fdr)*Ta;                                       %�����մ���?
nr = round(Tp*Fs);                      %�������ڵ���?
tr1 = 1:fix(nr)+1;                      
tr1 = tr1/Fs-Tp/2; 
Na = 8192;%2^14;                              %��λ�����
   
fclose(fd);

