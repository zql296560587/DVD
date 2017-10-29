%解析参数
fd = fopen('F:\P_20170609\PSAR_Stripmap_Br200MHz_Tr30us_Fsr260MHz_PRF1400_30K_20170609_1210_8_0.dat','rb');
% fd = fopen('F:\P_20170609\PSAR_Stripmap_Br30MHz_Tr30us_Fsr65MHz_PRF1400_16K_20170609_1246_18_0.dat','rb');
% fd = fopen('E:\P_20170524\CSAR_DP_X1V1_Br30MHz_Tr232us_Fsr51203MHz_20170524_1546_43_0.dat','rb');
% fd = fopen('I:\DWD\L波段数据\LSAR_Stripmap_Br200MHz_Tr20us_Fsr250MHz_PRF1000Hz_32K_20170609_1337_29_0.dat','rb');
% fseek(fd,(32768*4+128*8)*32768*2,'bof');

%% 解析包头
Ass = fread(fd,128,'uint8');
index = Ass(28);
Tmin = ( Ass(47) + Ass(48)*256 )*1e-6;      %采样延时
Tp = Ass(27,1)*1e-6;                        %脉冲宽度20us
% fc = 1300*1e6;             %L波段载频
fc = 390*1e6;                %P波段载频
Br = (Ass(34)*256+Ass(33))*1e6;              %带宽200M
Fs = (Ass(36)*256+Ass(35))*1e6;             %采样率250M
prf = (Ass(32)*256+Ass(31))/2;                  %PRF1000
Nr = (Ass(30)*256+Ass(29))*1024;            %距离向处理点数
H = (Ass(66)*256^3+Ass(65)*256^2+Ass(64)*256+Ass(63))*0.1;
% angle_elevation = Ass(74)*256+Ass(73);      %波束中心俯仰角
% if angle_elevation>0
%     angle_elevation = -1*angle_elevation*pi/180;
% end
% henggun = (Ass(50)*256+Ass(49))*0.001*pi/180;     %飞机横滚角
% fy = (Ass(52)*256+Ass(51))*0.001*pi/180;          %飞机俯仰角
% hangxiang = (Ass(54)*256+Ass(53))*0.01*pi/180;   %飞机航向角?
c = 3e8;                %光速
j = sqrt(-1);           %j^2=-1

Rmin = Tmin*c/2;
va = 400/3.6;           %速度  
% D = 1.36;             %L波段天线口径
D = 1.44;               %P波段天线口径

lamda = c/fc;           %波长
bin_r = c/2/Fs;         %距离向采样点间距离
gama = Br/Tp;            %距离向调频率?
delta_theta = 0.89*lamda/D;  %方位波束宽度
delta1 = rad2deg(delta_theta);

R = Rmin+((1:Nr)-1)*bin_r;      %各距离向作用距离


% angle_azimuth = pi/2;        %波束中心方位角
% [ReYzPit,ReYzAzim] = ReAirSpace(angle_azimuth, angle_elevation, henggun, fy, hangxiang);            %转到大地坐标系下
% thetas_c = acos(cos(ReYzAzim)*cos(ReYzPit));            %前斜角
% fdc = 2*va*cos(thetas_c)/lamda;                         %多普勒中心频率
% fdr = -2*va*va*sin(thetas_c)^2/lamda/Ra;                %多普勒调频率

Na = 16384;%25000;%2^14;                              %方位向点数
Kr=gama;                                %距离向调频率?
fr=(-Fs/2):Fs/(Nr):(Fs/2-Fs/(Nr));      %距离向频率
               
fclose(fd);
