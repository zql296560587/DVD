clc;
clear all ;
close all ;

SignOfChirpSlope = 1 ;  % chirp信号调频率的符号

FolderImageOutPut = 'E:\程序\jieguo\' ;        %输出路径
% FolderImageOutPut = 'E:\DWD程序\jieguoP1\' ;
         
if ~isdir(FolderImageOutPut)
    mkdir(FolderImageOutPut) ;
end

% setparameterP;            % P 波段解析惯导参数
setparameter;               % L 波段解析惯导参数

loop = 1; % L波段选取数据位置
%% 读取回波
% read_raw('F:\P_20170609\PSAR_Stripmap_Br200MHz_Tr30us_Fsr260MHz_PRF1400_30K_20170609_1210_8_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);
% read_raw('F:\P_20170609\PSAR_Stripmap_Br30MHz_Tr30us_Fsr65MHz_PRF1400_16K_20170609_1246_18_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);
read_raw('E:\L_20170609\LSAR_Stripmap_Br200MHz_Tr20us_Fsr250MHz_PRF1000Hz_32K_20170609_1205_4_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,16384,index,loop);         %读取当前回波数据
% read_raw('E:\P_20170524\CSAR_DP_X1V1_Br30MHz_Tr232us_Fsr51203MHz_20170524_1546_43_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);         %读取当前回波数据

%% 处理时修改的参数（为了和处理程序对应）
    % Nr = 15000;       % P 波段距离向处理选取点数
    Na = 8192;          %方位向处理点数
    %     Na = 25000/2-2;
    C = c ;         %光速
    Lambda = lamda;         %波长 
    Prf = prf;              %重频     
    NrNew = 8192;           %脉压后距离向点数（去除非完全积累点）

    FdcRefCalcu = 0 ;      % 雷达多普勒中心，目前是正侧
    Fdc_n = round(FdcRefCalcu/Prf);             %多普勒模糊数

    tr = Rmin*2/C + (0:Nr-1)'/Fs ;          % 距离快时间
    ta = (-Na/2:Na/2-1)'/Prf ;              % 方位慢时间

    fa = [0:Na/2-1 -Na/2:-1]'/Na*Prf ;      % 方位频率
    fr = [0:Nr/2-1 -Nr/2:-1]'/Nr*Fs ;       % 距离频率

    NrNewStart = 1;         %距离脉压后的数据起始位置
    % NrNewStart = 1200;    % P 波段（5月）距离脉压后的数据起始位置
    RminNew = Rmin+(NrNewStart-1)*(C/Fs/2) ;        %距离向去掉非完全积累点后的新的起始采样距离
    Rref = RminNew+NrNew/2*C/Fs/2 ;                 %参考距离（场景中心作用距离）
    RnNew = RminNew + (0:NrNew-1)'/Fs*C/2 ;         %距离向去掉非完全积累点后的新的作用距离
    frNew = [0:NrNew/2-1 -NrNew/2:-1]'/NrNew*Fs ;   %距离向去掉非完全积累点后的新的距离频率

    THETA1 = 0;     %前斜角，正侧，理论为0
    FdcCalcu = 0;   %fdc计算值，正侧，理论为0

    V = va;
    Beta = asind(H./RnNew) ; % 随斜距变化的擦地角
    FdrCalcu = -2*V^2*cosd(THETA1).^2/Lambda./RnNew;

%%
    disp('距离FFT，距离压缩，距离IFFT，去除距离向非完全积累点') ;
    StepRanComp ;
    ResultDisplayRaw(  'RanCompReal.dat' , 'RanCompImag.dat' , [FolderImageOutPut 'RanComp.raw'], Na , NrNew ) ;

%%
    disp('转置（由横向距离变为横向方位）') ;
    CornerTurn( 'RanCompReal.dat' , 'RanCompTurnReal.dat' , Na , NrNew , 4 ) ;
    CornerTurn( 'RanCompImag.dat' , 'RanCompTurnImag.dat' , Na , NrNew , 4 ) ;
    
%% 
    disp('方位向去高频分量') ;
    StepAzi_high_delete ;
    
%% 
    disp('杂波锁定，估计fdc');
    StepClutterLock ;           %估计多普勒中心频率fdc
    
%%
    disp('转置（由横向方位变为横向距离）') ;
    CornerTurn( 'RanCompTurnReal1.dat' , 'RanCompReal1.dat' , NrNew , Na , 4 ) ;
    CornerTurn( 'RanCompTurnImag1.dat' , 'RanCompImag1.dat' , NrNew , Na , 4 ) ;
    
%%
    disp('残余走动校正，去多普勒中心');
    StepResiRanWalkCorr;
    ResultDisplayRaw( 'ResiRanWalkCorrReal.dat' , 'ResiRanWalkCorrImag.dat' ,  'StepResiRanWalkCorr.raw', Na , NrNew ) ;

%%
    disp('转置（由横向距离转为横向方位）') ;
    CornerTurn( 'ResiRanWalkCorrReal.dat' , 'ResiRanWalkCorrTurnReal.dat' , Na , NrNew , 4 ) ;
    CornerTurn( 'ResiRanWalkCorrImag.dat' , 'ResiRanWalkCorrTurnImag.dat' , Na , NrNew , 4 ) ;

%%
    disp('方位FFT') ;
    StepAziFft ;
    ResultDisplayRaw( 'AziFftReal.dat' , 'AziFftImag.dat' , 'AziFFT.raw', NrNew , Na ) ;

%%
    disp('转置（由横向方位转为横向距离）') ;
    CornerTurn( 'AziFftReal.dat' , 'AziFftTurnReal.dat' , NrNew , Na , 2 ) ;
    CornerTurn( 'AziFftImag.dat' , 'AziFftTurnImag.dat' , NrNew , Na , 2 ) ;
    ResultDisplayRaw( 'AziFftTurnReal.dat' , 'AziFftTurnImag.dat' ,  'AziFftTurn.raw', Na , NrNew ) ;

%%
    disp('距离FFT，距离弯曲校正，距离IFFT') ;
    StepRanCurveCorr ;
    ResultDisplayRaw( 'RanCurveCorrReal.dat' , 'RanCurveCorrImag.dat' ,  'RanCurveCorr.raw', Na , NrNew ) ;

%%
    disp('转置（由横向距离转为横向方位）') ;
    CornerTurn( 'RanCurveCorrReal.dat' , 'RanCurveCorrTurnReal.dat' , Na , NrNew , 2 ) ;
    CornerTurn( 'RanCurveCorrImag.dat' , 'RanCurveCorrTurnImag.dat' , Na , NrNew , 2 ) ;

%%
    disp('方位IFFT') ;
    StepAziIFft ;

%%
    disp('自聚焦估计fdr') ;
    StepAutoFocus ;         %估计多普勒调频率fdr

%%
    disp('方位向处理') ;
    StepAziComp_fenkuai;
    
%%  结果
    ResultDisplayRaw( 'AziCompReal.dat' , 'AziCompImag.dat' ,  [FolderImageOutPut 'Azix3Comp'  '_' num2str(Na) '_' num2str(loop+3) '_.raw'], NrNew , N_Ba ) ;

