%% �Ծ۽�

disp('�Ծ۽������Աȶȷ����ƶ����յ�Ƶ��') ;

% FdrEstimDepth = 4 ;  % ÿFdrEstimDepth�������Ź��Ƴ�һ��Fdrֵ��Ҫ���ܱ������Ÿ���NrNew����
FdrEstimDepth = 1 ;  % ÿFdrEstimDepth�������Ź��Ƴ�һ��Fdrֵ��Ҫ���ܱ������Ÿ���NrNew����
FdrStepInit = 50 ;      % Fdr��ʼ��������
FdrStepFinal = 1 ;      % Fdr������������
MaxIterateNum = 40 ;    % ����������


% DiffCurvePathCompenPhase = CurvePathCompenPhase(2:end) - CurvePathCompenPhase(1:end-1) ;
% DiffDiffCurvePathCompenPhase = DiffCurvePathCompenPhase(2:end) - DiffCurvePathCompenPhase(1:end-1) ;
% SecondCoeffCurveCompen = mean(DiffDiffCurvePathCompenPhase)*Prf*Prf/2/pi ;  % ����������λ�ж�����λ��ϵ��


% FdrInit = -2*V^2*sin(ThetaVEstim(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)).^2/Lambda ...
%     ./(RnNew(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)) ;  % ��ʼFdr
% FdrInit = -2*V^2*sind(ThetaVEstim(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)).^2/Lambda ...
%     ./(RnNew(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)) ;  % �����Ǽ��ٶȵĳ�ʼFdr

% FdrInit = -2*V^2*sind(ThetaVEstim(NrNStart:FdrEstimDepth:NrNEnd)).^2/Lambda ...
%     ./(RnNew(NrNStart:FdrEstimDepth:NrNEnd)) ;  % �����Ǽ��ٶȵĳ�ʼFdr

% FdrInit = FdrCalcu(NrNStart:FdrEstimDepth:NrNEnd);

FdrInit = -2*V^2*(1-(Lambda*FdcEstim(1:FdrEstimDepth:NrNew)/2/V).^2)./Lambda./(RnNew(1:FdrEstimDepth:NrNew));    % �����յ�Ƶ����б��ı仯

FdrMaxContrast ;

%% 
disp('�Ծ۽����Թ��ƽ�������˲������') ;

% ��RΪ��������һ�����
M = 1:FdrEstimDepth:NrNew;
X = RnNew(1:FdrEstimDepth:NrNew) ;
Y = FdrOfMaxContrast;
CoeffOfFdr = polyfit( X , Y , 1 ) ;
FdrEstim1 =  CoeffOfFdr(1).*X + CoeffOfFdr(2)  ;
% ȥ������3sigma��ֵ���������������
RemoveIndex = find(abs(FdrOfMaxContrast-FdrEstim1) ...
    <std(FdrOfMaxContrast-FdrEstim1)) ;
X = X(RemoveIndex) ;
Y = Y(RemoveIndex) ;
CoeffOfFdr = polyfit( X , Y , 2 ) ;
FdrEstim =  CoeffOfFdr(1).*RnNew.^2 + CoeffOfFdr(2).*RnNew + CoeffOfFdr(3) ;


figure; plot( (1:FdrEstimDepth:NrNew)' , FdrInit ,'b' ) ;
hold on , plot( (1:FdrEstimDepth:NrNew)' , FdrOfMaxContrast , 'r' ) ;
plot(FdrCalcu , 'k' ) ;
plot( FdrEstim , 'g' ) ;
xlabel('������') , ylabel('�����յ�Ƶ�ʣ�Hz/s��') ;
title('�����յ�Ƶ�ʳ�ʼֵ,����ֵ,���ֵ') ;
legend( 'Fdr��ʼֵ' , 'Fdr����ֵ' , 'Fdr�ߵ�ֵ' , 'Fdr���ֵ' , 'location' , 'best' ) ;

figure ; plot( FdrCalcu , 'r' ) ;
hold on  ;
plot( FdrEstim , 'g' ) ;
xlabel('������') , ylabel('�����յ�Ƶ�ʣ�Hz/s��') ;
title('�����յ�Ƶ�ʹߵ�ֵ,���ֵ') ;
legend( 'Fdr�ߵ�ֵ' , 'Fdr���ֵ' , 'location' , 'best' ) ;

