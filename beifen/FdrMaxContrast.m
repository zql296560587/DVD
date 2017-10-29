%% ���Աȶȷ����ƶ����յ�Ƶ��

FdrOfMaxContrast = zeros(NrNew/FdrEstimDepth,1) ;    % ���Աȶȷ��õ���Fdr����ֵ
FirstAutofocusSuccess = -1 ;    % ��һ���Ծ۽��ɹ��ľ۽����
IsAutofocusSuccess = true(NrNew/FdrEstimDepth,1) ;   % ÿ���۽�����Ƿ��Ծ۽��ɹ��ı�־
% WinAzi = hamming(Na) ;
WinAzi = ones(Na,1) ;

FidReadReal = fopen( 'AziIFftReal.dat' , 'r' ) ;
FidReadImag = fopen( 'AziIFftImag.dat' , 'r' ) ;

for m = 1 : NrNew/FdrEstimDepth
    
%     m
    
    %% 1��������ÿ���۽�������ҵ��Աȶ����ľ�����
    
    FdcOfTheDepth = FdcEstim(m) - FdcRefEstim ;
    MaxContrast = 0 ;   % ��¼���ĶԱȶ�
    Data = fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' ) ;
    DataOfMaxContrast = (Data);
     
    %% 2��Ȼ�����Աȶ��Ծ۽�
    
    FdrCurrent = FdrInit(m) ;       % ��ǰ��Fdr
    IterateNumCurrent = 0 ;         % ��������
    FdrStepCurrent = FdrStepInit ;  % ��ǰ�ĵ�������
    ThreeDataContrast = zeros(3,1) ;    % ����Fdr��ѹ�õ��������Աȶ�
    
    % ��1�������ʼFdr�Լ��Ӽ�һ����ʼ���������ĶԱȶ�
    IterateNumCurrent = IterateNumCurrent + 1 ;
    ThreeDataContrast(2) = MaxContrast ;
    for k = -1 : 2 : 1
%         Hazi = exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
        Hazi =  WinAzi.*exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
        DataComp = abs( fft( DataOfMaxContrast .* Hazi ) ) ;
        ThreeDataContrast(k+2) = sqrt(var(DataComp))/mean(DataComp);        
    end
    [NoUse, IndexOfMaxContrast] = max(ThreeDataContrast) ;
    
    % ��2�����һֱ�ǵ�һ������һֱ�ǵ�����Fdr�ĶԱȶ�����򲽳�����
    while IndexOfMaxContrast~=2
        
        IterateNumCurrent = IterateNumCurrent + 1;
        
        if IterateNumCurrent>MaxIterateNum
            % ������������������������������Ծ۽�ʧ��
            FdrOfMaxContrast(m) = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent ;
            IsAutofocusSuccess(m) = false ;
            break ;
        end
        
        FdrCurrent = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent;
        ThreeDataContrast(2) = ThreeDataContrast(IndexOfMaxContrast) ;
        
        for k = (IndexOfMaxContrast-2) : 1 : (IndexOfMaxContrast-2)
%             Hazi = exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
            Hazi =  WinAzi.*exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
            DataComp = abs( fft( DataOfMaxContrast .* Hazi ) ) ;
            ThreeDataContrast(k+2) = sqrt(var(DataComp))/mean(DataComp);
        end
        [NoUse, IndexOfMaxContrast] = max(ThreeDataContrast) ;
        
    end
    
    % ��3��һ�����ֵڶ���Fdr�ĶԱȶ�����򲽳�����
    while true
        
        IterateNumCurrent = IterateNumCurrent + 1;

        if  IterateNumCurrent>MaxIterateNum
            % ������������������������������Ծ۽�ʧ��
            FdrOfMaxContrast(m) = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent ;
            IsAutofocusSuccess(m) = false ;
            break ;
        end
        
        if FdrStepCurrent<=FdrStepFinal
            % �����ǰ����С�����ղ��������Ծ۽��ɹ�����¼��ǰFdr
            FdrOfMaxContrast(m) = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent ;
            if FirstAutofocusSuccess==-1
                % ������ǵ�һ���Ծ۽��ɹ��ľ۽���ȣ����¼����
                FirstAutofocusSuccess = m ;
            end
            break ;
        end
        
        if IndexOfMaxContrast~=2
            % �����һ���������Fdr�ĶԱȶ��������µ�ǰFdr
            FdrCurrent = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent;
            ThreeDataContrast(2) = ThreeDataContrast(IndexOfMaxContrast);
        end
        
        FdrStepCurrent = FdrStepCurrent * 0.5;
        
        for k = -1 : 2 : 1
%             Hazi = exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
            Hazi =  WinAzi.*exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
            DataComp = abs( fft( DataOfMaxContrast .* Hazi ) ) ;
            ThreeDataContrast(k+2) = sqrt(var(DataComp))/mean(DataComp);
        end
        [NoUse, IndexOfMaxContrast] = max(ThreeDataContrast) ;
        
    end
        
end

fclose all ;

%% 3�����Ծ۽�ʧ�ܵľ۽���Ƚ��д���
if(FirstAutofocusSuccess == -1)
    % ���ȫ���Ծ۽�ʧ�ܣ���ʹ�ó�ʼֵ
    FdrOfMaxContrast = FdrInit ;
else
    % ���û��ȫ���Ծ۽�ʧ��
    % ���ڵ�һ���۽���ȣ�ʹ�õ�һ���Ծ۽��ɹ���Fdrֵ
    FdrOfMaxContrast(1) = FdrOfMaxContrast(FirstAutofocusSuccess);
    % ���ڵڶ����������ڶ����۽�������Ծ۽�ʧ�ܵľ۽���ȣ������һ�۽����
    % �Ծ۽��ɹ���ȡǰһ�����һ���ľ�ֵ������ȡǰһ��ֵ��
    for k = 2 : NrNew/FdrEstimDepth-1
        if ~IsAutofocusSuccess(k)
            if IsAutofocusSuccess(k+1)
                FdrOfMaxContrast(k) = (FdrOfMaxContrast(k-1)+FdrOfMaxContrast(k+1))/2;
            else
                FdrOfMaxContrast(k) = FdrOfMaxContrast(k-1);
            end
        end
    end
    % �������һ���۽���ȣ�����Ծ۽�ʧ�ܣ���ȡǰһ��ֵ
    k = k+1 ;
    if ~IsAutofocusSuccess(k)
        FdrOfMaxContrast(k) = FdrOfMaxContrast(k-1);
    end    
end

