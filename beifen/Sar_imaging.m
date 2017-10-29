clc;
clear all ;
close all ;

SignOfChirpSlope = 1 ;  % chirp�źŵ�Ƶ�ʵķ���

FolderImageOutPut = 'E:\����\jieguo\' ;        %���·��
% FolderImageOutPut = 'E:\DWD����\jieguoP1\' ;
         
if ~isdir(FolderImageOutPut)
    mkdir(FolderImageOutPut) ;
end

% setparameterP;            % P ���ν����ߵ�����
setparameter;               % L ���ν����ߵ�����

loop = 1; % L����ѡȡ����λ��
%% ��ȡ�ز�
% read_raw('F:\P_20170609\PSAR_Stripmap_Br200MHz_Tr30us_Fsr260MHz_PRF1400_30K_20170609_1210_8_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);
% read_raw('F:\P_20170609\PSAR_Stripmap_Br30MHz_Tr30us_Fsr65MHz_PRF1400_16K_20170609_1246_18_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);
read_raw('E:\L_20170609\LSAR_Stripmap_Br200MHz_Tr20us_Fsr250MHz_PRF1000Hz_32K_20170609_1205_4_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,16384,index,loop);         %��ȡ��ǰ�ز�����
% read_raw('E:\P_20170524\CSAR_DP_X1V1_Br30MHz_Tr232us_Fsr51203MHz_20170524_1546_43_0.dat',[FolderImageOutPut 'EchoReal.dat'],[FolderImageOutPut 'EchoImag.dat'],Na,Nr,index);         %��ȡ��ǰ�ز�����

%% ����ʱ�޸ĵĲ�����Ϊ�˺ʹ�������Ӧ��
    % Nr = 15000;       % P ���ξ�������ѡȡ����
    Na = 8192;          %��λ�������
    %     Na = 25000/2-2;
    C = c ;         %����
    Lambda = lamda;         %���� 
    Prf = prf;              %��Ƶ     
    NrNew = 8192;           %��ѹ������������ȥ������ȫ���۵㣩

    FdcRefCalcu = 0 ;      % �״���������ģ�Ŀǰ������
    Fdc_n = round(FdcRefCalcu/Prf);             %������ģ����

    tr = Rmin*2/C + (0:Nr-1)'/Fs ;          % �����ʱ��
    ta = (-Na/2:Na/2-1)'/Prf ;              % ��λ��ʱ��

    fa = [0:Na/2-1 -Na/2:-1]'/Na*Prf ;      % ��λƵ��
    fr = [0:Nr/2-1 -Nr/2:-1]'/Nr*Fs ;       % ����Ƶ��

    NrNewStart = 1;         %������ѹ���������ʼλ��
    % NrNewStart = 1200;    % P ���Σ�5�£�������ѹ���������ʼλ��
    RminNew = Rmin+(NrNewStart-1)*(C/Fs/2) ;        %������ȥ������ȫ���۵����µ���ʼ��������
    Rref = RminNew+NrNew/2*C/Fs/2 ;                 %�ο����루�����������þ��룩
    RnNew = RminNew + (0:NrNew-1)'/Fs*C/2 ;         %������ȥ������ȫ���۵����µ����þ���
    frNew = [0:NrNew/2-1 -NrNew/2:-1]'/NrNew*Fs ;   %������ȥ������ȫ���۵����µľ���Ƶ��

    THETA1 = 0;     %ǰб�ǣ����࣬����Ϊ0
    FdcCalcu = 0;   %fdc����ֵ�����࣬����Ϊ0

    V = va;
    Beta = asind(H./RnNew) ; % ��б��仯�Ĳ��ؽ�
    FdrCalcu = -2*V^2*cosd(THETA1).^2/Lambda./RnNew;

%%
    disp('����FFT������ѹ��������IFFT��ȥ�����������ȫ���۵�') ;
    StepRanComp ;
    ResultDisplayRaw(  'RanCompReal.dat' , 'RanCompImag.dat' , [FolderImageOutPut 'RanComp.raw'], Na , NrNew ) ;

%%
    disp('ת�ã��ɺ�������Ϊ����λ��') ;
    CornerTurn( 'RanCompReal.dat' , 'RanCompTurnReal.dat' , Na , NrNew , 4 ) ;
    CornerTurn( 'RanCompImag.dat' , 'RanCompTurnImag.dat' , Na , NrNew , 4 ) ;
    
%% 
    disp('��λ��ȥ��Ƶ����') ;
    StepAzi_high_delete ;
    
%% 
    disp('�Ӳ�����������fdc');
    StepClutterLock ;           %���ƶ���������Ƶ��fdc
    
%%
    disp('ת�ã��ɺ���λ��Ϊ������룩') ;
    CornerTurn( 'RanCompTurnReal1.dat' , 'RanCompReal1.dat' , NrNew , Na , 4 ) ;
    CornerTurn( 'RanCompTurnImag1.dat' , 'RanCompImag1.dat' , NrNew , Na , 4 ) ;
    
%%
    disp('�����߶�У����ȥ����������');
    StepResiRanWalkCorr;
    ResultDisplayRaw( 'ResiRanWalkCorrReal.dat' , 'ResiRanWalkCorrImag.dat' ,  'StepResiRanWalkCorr.raw', Na , NrNew ) ;

%%
    disp('ת�ã��ɺ������תΪ����λ��') ;
    CornerTurn( 'ResiRanWalkCorrReal.dat' , 'ResiRanWalkCorrTurnReal.dat' , Na , NrNew , 4 ) ;
    CornerTurn( 'ResiRanWalkCorrImag.dat' , 'ResiRanWalkCorrTurnImag.dat' , Na , NrNew , 4 ) ;

%%
    disp('��λFFT') ;
    StepAziFft ;
    ResultDisplayRaw( 'AziFftReal.dat' , 'AziFftImag.dat' , 'AziFFT.raw', NrNew , Na ) ;

%%
    disp('ת�ã��ɺ���λתΪ������룩') ;
    CornerTurn( 'AziFftReal.dat' , 'AziFftTurnReal.dat' , NrNew , Na , 2 ) ;
    CornerTurn( 'AziFftImag.dat' , 'AziFftTurnImag.dat' , NrNew , Na , 2 ) ;
    ResultDisplayRaw( 'AziFftTurnReal.dat' , 'AziFftTurnImag.dat' ,  'AziFftTurn.raw', Na , NrNew ) ;

%%
    disp('����FFT����������У��������IFFT') ;
    StepRanCurveCorr ;
    ResultDisplayRaw( 'RanCurveCorrReal.dat' , 'RanCurveCorrImag.dat' ,  'RanCurveCorr.raw', Na , NrNew ) ;

%%
    disp('ת�ã��ɺ������תΪ����λ��') ;
    CornerTurn( 'RanCurveCorrReal.dat' , 'RanCurveCorrTurnReal.dat' , Na , NrNew , 2 ) ;
    CornerTurn( 'RanCurveCorrImag.dat' , 'RanCurveCorrTurnImag.dat' , Na , NrNew , 2 ) ;

%%
    disp('��λIFFT') ;
    StepAziIFft ;

%%
    disp('�Ծ۽�����fdr') ;
    StepAutoFocus ;         %���ƶ����յ�Ƶ��fdr

%%
    disp('��λ����') ;
    StepAziComp_fenkuai;
    
%%  ���
    ResultDisplayRaw( 'AziCompReal.dat' , 'AziCompImag.dat' ,  [FolderImageOutPut 'Azix3Comp'  '_' num2str(Na) '_' num2str(loop+3) '_.raw'], NrNew , N_Ba ) ;

