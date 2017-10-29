%% �Ӳ�����

% FdcRefCalcu�����ݹߵ�����Ĳο�б�ദ�Ķ��������ģ�
% FdcOfTimeCorr��ʱ����ط����Ƴ�����б��仯�ļ�ȥ��FdcRefCalcu�Ķ��������ģ�
% FdcOfCorr��ʱ����ط����Ƴ�����б��仯�ļ�ȥ��FdcRefCalcu�Ķ��������Ľ���ƽ����
% FdcEstim��FdcOfCorr����FdcRefCalcu��������͵Ľ���������չ��Ƶõ��Ķ��������ģ�
% FdcRefEstim��FdcEstim�ο�б���ֵ��

%% ʱ����ط����ƶ���������
disp('�Ӳ�������ʱ����ط����ƶ���������') ;

FdcEstimDepth = 1 ;    % ÿFdcEstimDepth�������Ź���һ�Σ�Ҫ���ܱ������Ÿ���NrNew����
FdcTimeCorrelate ;

%% 3���Թ��ƽ�������
disp('�Ӳ��������Թ��ƽ�����н����');

LenOfFdc = length(FdcOfTimeCorr) ;
SigFdc = exp(1j*2*pi*FdcOfTimeCorr/Prf) ;
% Nfft2 = LenOfFdc*10 ;
Nfft2 = Nr ;

Temp2 = fft( conj(SigFdc(1:LenOfFdc/2)).*SigFdc((LenOfFdc/2+1):LenOfFdc) , Nfft2 ) ;
[Unused,MaxIndex2] = max(abs(Temp2)) ;
if MaxIndex2>Nfft2/2
    MaxIndex2 = MaxIndex2-Nfft2 ;
end
FdcTemp2 = (MaxIndex2-1)*1/Nfft2/(LenOfFdc/2)*Prf/2 ;           %�õ�fdc�Ķ���ϵ��

Nfft1 = Nr ;
Temp1 = fft( SigFdc.*exp(-2j*pi*FdcTemp2*(0:LenOfFdc-1)'.^2/Prf) , Nfft1 ) ;
[Unused,MaxIndex1] = max(abs(Temp1)) ;
MaxValue1 = Temp1(MaxIndex1) ;
if MaxIndex1>Nfft1/2
    MaxIndex1 = MaxIndex1-Nfft1 ;
end
FdcTemp1 = (MaxIndex1-1)*1/Nfft1*Prf ;          %�õ�fdc��һ��ϵ��
FdcTemp0 = angle(MaxValue1)/2/pi*Prf ;          %�õ�fdc�ĳ�����
FdcTemp = FdcTemp0+FdcTemp1*(0:LenOfFdc-1)'+FdcTemp2*(0:LenOfFdc-1)'.^2 ;       

FdcOfCorr = FdcOfTimeCorr + round((FdcTemp-FdcOfTimeCorr)/Prf)*Prf ;

% ������ģ����������
FdcDeAmb = FdcOfCorr + Fdc_n*Prf ;
% FdcDeAmb = FdcOfCorr + FdcRefCalcu ;
FdcDeAmb = FdcDeAmb - Prf*round( (mean(FdcDeAmb)-FdcRefCalcu)/Prf ) ;           %����ƽ��

%% 4���Խ���ƽ�������˲�����ϴ���
disp('�Ӳ��������Խ���ƽ�������˲������') ;

% % % ��RΪ��������һ�����
% X = RnNew;
% Y = FdcDeAmb ;
% CoeffOfFdc = polyfit( X , Y , 1 ) ;
% FdcEstim = CoeffOfFdc(1).*X+CoeffOfFdc(2) ;
% 
% PolyfitError = FdcEstim - FdcDeAmb ;
% SigmaPolyfitError = std(PolyfitError) ;
% IndexPolyfit = abs(PolyfitError)<1*SigmaPolyfitError ;  % �˲����������
% X = X(IndexPolyfit) ;
% Y = Y(IndexPolyfit) ;
% 
% CoeffOfFdc = polyfit( X , Y , 2 ) ;
% FdcEstim = CoeffOfFdc(1).*RnNew.^2+CoeffOfFdc(2).*RnNew+CoeffOfFdc(3) ;  

FdcEstim = ones(NrNew,1)*mean(FdcDeAmb);        %���ڶ�Ŀ��������о�������ˣ�fdcѡȡ��Ŀ�������������ľ��봦��һ����ֵ
 
SinTheta = Lambda*FdcEstim/2/V;

if abs(SinTheta(1)) >1
    SinTheta = ones(1,NrNew);
end
CosTheta = sqrt(1-SinTheta.^2);
ThetaVEstim = acosd(CosTheta) ;
FdcRefEstim= FdcEstim(NrNew/2) ;
ThetaVRefEstim = asind(Lambda*FdcRefEstim/2/V) ;         % ���ù��ƵĶ���������ֻ��ǰб�ǽ����˸���

%%
figure;plot((1:FdcEstimDepth:NrNew)' , FdcOfTimeCorr);
hold on , plot((1:FdcEstimDepth:NrNew)' , FdcOfCorr , 'r');
xlabel('������') , ylabel('���������ģ�Hz��') ; axis tight ;
title( '��б��仯�Ķ���������' ) , legend('Fdc����ֵ','Fdc�����ֵ','location','best') ;

figure; plot( (1:FdcEstimDepth:NrNew)' , FdcDeAmb ) ;
hold on , plot( (1:NrNew)' , FdcEstim , 'r' ) ;
plot( (1:NrNew)' , FdcCalcu ,'g' ) ;
xlabel('������') , ylabel('���������ģ�Hz��') ; axis tight ;
title( '��б��仯�Ķ���������' ) , legend('Fdc����ֵ','Fdc���ֵ','Fdc�ߵ�ֵ','location','best') ;
