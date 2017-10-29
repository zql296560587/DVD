%% 杂波锁定

% FdcRefCalcu：根据惯导计算的参考斜距处的多普勒中心；
% FdcOfTimeCorr：时域相关法估计出的随斜距变化的减去了FdcRefCalcu的多普勒中心；
% FdcOfCorr：时域相关法估计出的随斜距变化的减去了FdcRefCalcu的多普勒中心解缠绕结果；
% FdcEstim：FdcOfCorr加上FdcRefCalcu并曲线拟和的结果，即最终估计得到的多普勒中心；
% FdcRefEstim：FdcEstim参考斜距的值；

%% 时域相关法估计多普勒中心
disp('杂波锁定：时域相关法估计多普勒中心') ;

FdcEstimDepth = 1 ;    % 每FdcEstimDepth个距离门估计一次，要求能被距离门个数NrNew整除
FdcTimeCorrelate ;

%% 3、对估计结果解缠绕
disp('杂波锁定：对估计结果进行解缠绕');

LenOfFdc = length(FdcOfTimeCorr) ;
SigFdc = exp(1j*2*pi*FdcOfTimeCorr/Prf) ;
% Nfft2 = LenOfFdc*10 ;
Nfft2 = Nr ;

Temp2 = fft( conj(SigFdc(1:LenOfFdc/2)).*SigFdc((LenOfFdc/2+1):LenOfFdc) , Nfft2 ) ;
[Unused,MaxIndex2] = max(abs(Temp2)) ;
if MaxIndex2>Nfft2/2
    MaxIndex2 = MaxIndex2-Nfft2 ;
end
FdcTemp2 = (MaxIndex2-1)*1/Nfft2/(LenOfFdc/2)*Prf/2 ;           %得到fdc的二次系数

Nfft1 = Nr ;
Temp1 = fft( SigFdc.*exp(-2j*pi*FdcTemp2*(0:LenOfFdc-1)'.^2/Prf) , Nfft1 ) ;
[Unused,MaxIndex1] = max(abs(Temp1)) ;
MaxValue1 = Temp1(MaxIndex1) ;
if MaxIndex1>Nfft1/2
    MaxIndex1 = MaxIndex1-Nfft1 ;
end
FdcTemp1 = (MaxIndex1-1)*1/Nfft1*Prf ;          %得到fdc的一次系数
FdcTemp0 = angle(MaxValue1)/2/pi*Prf ;          %得到fdc的常数项
FdcTemp = FdcTemp0+FdcTemp1*(0:LenOfFdc-1)'+FdcTemp2*(0:LenOfFdc-1)'.^2 ;       

FdcOfCorr = FdcOfTimeCorr + round((FdcTemp-FdcOfTimeCorr)/Prf)*Prf ;

% 多普勒模糊数的修正
FdcDeAmb = FdcOfCorr + Fdc_n*Prf ;
% FdcDeAmb = FdcOfCorr + FdcRefCalcu ;
FdcDeAmb = FdcDeAmb - Prf*round( (mean(FdcDeAmb)-FdcRefCalcu)/Prf ) ;           %解缠绕结果

%% 4、对解缠绕结果进行滤波和拟合处理
disp('杂波锁定：对解缠绕结果进行滤波和拟合') ;

% % % 以R为变量进行一次拟合
% X = RnNew;
% Y = FdcDeAmb ;
% CoeffOfFdc = polyfit( X , Y , 1 ) ;
% FdcEstim = CoeffOfFdc(1).*X+CoeffOfFdc(2) ;
% 
% PolyfitError = FdcEstim - FdcDeAmb ;
% SigmaPolyfitError = std(PolyfitError) ;
% IndexPolyfit = abs(PolyfitError)<1*SigmaPolyfitError ;  % 滤波后重新拟合
% X = X(IndexPolyfit) ;
% Y = Y(IndexPolyfit) ;
% 
% CoeffOfFdc = polyfit( X , Y , 2 ) ;
% FdcEstim = CoeffOfFdc(1).*RnNew.^2+CoeffOfFdc(2).*RnNew+CoeffOfFdc(3) ;  

FdcEstim = ones(NrNew,1)*mean(FdcDeAmb);        %由于对目标区域进行精成像，因此，fdc选取与目标所在区域（中心距离处）一样的值
 
SinTheta = Lambda*FdcEstim/2/V;

if abs(SinTheta(1)) >1
    SinTheta = ones(1,NrNew);
end
CosTheta = sqrt(1-SinTheta.^2);
ThetaVEstim = acosd(CosTheta) ;
FdcRefEstim= FdcEstim(NrNew/2) ;
ThetaVRefEstim = asind(Lambda*FdcRefEstim/2/V) ;         % 利用估计的多普勒中心只对前斜角进行了更新

%%
figure;plot((1:FdcEstimDepth:NrNew)' , FdcOfTimeCorr);
hold on , plot((1:FdcEstimDepth:NrNew)' , FdcOfCorr , 'r');
xlabel('距离门') , ylabel('多普勒中心（Hz）') ; axis tight ;
title( '随斜距变化的多普勒中心' ) , legend('Fdc估计值','Fdc解缠绕值','location','best') ;

figure; plot( (1:FdcEstimDepth:NrNew)' , FdcDeAmb ) ;
hold on , plot( (1:NrNew)' , FdcEstim , 'r' ) ;
plot( (1:NrNew)' , FdcCalcu ,'g' ) ;
xlabel('距离门') , ylabel('多普勒中心（Hz）') ; axis tight ;
title( '随斜距变化的多普勒中心' ) , legend('Fdc估计值','Fdc拟合值','Fdc惯导值','location','best') ;
