%% 自聚焦

disp('自聚焦：最大对比度法估计多普勒调频率') ;

% FdrEstimDepth = 4 ;  % 每FdrEstimDepth个距离门估计出一个Fdr值，要求能被距离门个数NrNew整除
FdrEstimDepth = 1 ;  % 每FdrEstimDepth个距离门估计出一个Fdr值，要求能被距离门个数NrNew整除
FdrStepInit = 50 ;      % Fdr初始搜索步长
FdrStepFinal = 1 ;      % Fdr最终搜索步长
MaxIterateNum = 40 ;    % 最大迭代次数


% DiffCurvePathCompenPhase = CurvePathCompenPhase(2:end) - CurvePathCompenPhase(1:end-1) ;
% DiffDiffCurvePathCompenPhase = DiffCurvePathCompenPhase(2:end) - DiffCurvePathCompenPhase(1:end-1) ;
% SecondCoeffCurveCompen = mean(DiffDiffCurvePathCompenPhase)*Prf*Prf/2/pi ;  % 弯曲补偿相位中二次相位的系数


% FdrInit = -2*V^2*sin(ThetaVEstim(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)).^2/Lambda ...
%     ./(RnNew(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)) ;  % 初始Fdr
% FdrInit = -2*V^2*sind(ThetaVEstim(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)).^2/Lambda ...
%     ./(RnNew(floor(FdrEstimDepth/2)+1:FdrEstimDepth:NrNew)) ;  % 不考虑加速度的初始Fdr

% FdrInit = -2*V^2*sind(ThetaVEstim(NrNStart:FdrEstimDepth:NrNEnd)).^2/Lambda ...
%     ./(RnNew(NrNStart:FdrEstimDepth:NrNEnd)) ;  % 不考虑加速度的初始Fdr

% FdrInit = FdrCalcu(NrNStart:FdrEstimDepth:NrNEnd);

FdrInit = -2*V^2*(1-(Lambda*FdcEstim(1:FdrEstimDepth:NrNew)/2/V).^2)./Lambda./(RnNew(1:FdrEstimDepth:NrNew));    % 多普勒调频率随斜距的变化

FdrMaxContrast ;

%% 
disp('自聚焦：对估计结果进行滤波和拟合') ;

% 以R为变量进行一次拟合
M = 1:FdrEstimDepth:NrNew;
X = RnNew(1:FdrEstimDepth:NrNew) ;
Y = FdrOfMaxContrast;
CoeffOfFdr = polyfit( X , Y , 1 ) ;
FdrEstim1 =  CoeffOfFdr(1).*X + CoeffOfFdr(2)  ;
% 去掉超过3sigma的值，并重新曲线拟和
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
xlabel('距离门') , ylabel('多普勒调频率（Hz/s）') ;
title('多普勒调频率初始值,估计值,拟合值') ;
legend( 'Fdr初始值' , 'Fdr估计值' , 'Fdr惯导值' , 'Fdr拟合值' , 'location' , 'best' ) ;

figure ; plot( FdrCalcu , 'r' ) ;
hold on  ;
plot( FdrEstim , 'g' ) ;
xlabel('距离门') , ylabel('多普勒调频率（Hz/s）') ;
title('多普勒调频率惯导值,拟合值') ;
legend( 'Fdr惯导值' , 'Fdr拟合值' , 'location' , 'best' ) ;

