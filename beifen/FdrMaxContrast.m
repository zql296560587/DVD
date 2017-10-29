%% 最大对比度法估计多普勒调频率

FdrOfMaxContrast = zeros(NrNew/FdrEstimDepth,1) ;    % 最大对比度法得到的Fdr估计值
FirstAutofocusSuccess = -1 ;    % 第一个自聚焦成功的聚焦深度
IsAutofocusSuccess = true(NrNew/FdrEstimDepth,1) ;   % 每个聚焦深度是否自聚焦成功的标志
% WinAzi = hamming(Na) ;
WinAzi = ones(Na,1) ;

FidReadReal = fopen( 'AziIFftReal.dat' , 'r' ) ;
FidReadImag = fopen( 'AziIFftImag.dat' , 'r' ) ;

for m = 1 : NrNew/FdrEstimDepth
    
%     m
    
    %% 1、首先在每个聚焦深度内找到对比度最大的距离门
    
    FdcOfTheDepth = FdcEstim(m) - FdcRefEstim ;
    MaxContrast = 0 ;   % 记录最大的对比度
    Data = fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' ) ;
    DataOfMaxContrast = (Data);
     
    %% 2、然后最大对比度自聚焦
    
    FdrCurrent = FdrInit(m) ;       % 当前的Fdr
    IterateNumCurrent = 0 ;         % 迭代次数
    FdrStepCurrent = FdrStepInit ;  % 当前的迭代步长
    ThreeDataContrast = zeros(3,1) ;    % 三个Fdr脉压得到的三个对比度
    
    % （1）计算初始Fdr以及加减一个初始搜索步长的对比度
    IterateNumCurrent = IterateNumCurrent + 1 ;
    ThreeDataContrast(2) = MaxContrast ;
    for k = -1 : 2 : 1
%         Hazi = exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
        Hazi =  WinAzi.*exp( -2j*pi*FdcOfTheDepth*ta -1j*pi*(FdrCurrent+k*FdrStepCurrent)*ta.^2 ) ;
        DataComp = abs( fft( DataOfMaxContrast .* Hazi ) ) ;
        ThreeDataContrast(k+2) = sqrt(var(DataComp))/mean(DataComp);        
    end
    [NoUse, IndexOfMaxContrast] = max(ThreeDataContrast) ;
    
    % （2）如果一直是第一个或者一直是第三个Fdr的对比度最大，则步长不变
    while IndexOfMaxContrast~=2
        
        IterateNumCurrent = IterateNumCurrent + 1;
        
        if IterateNumCurrent>MaxIterateNum
            % 如果迭代次数超过最大迭代次数，则自聚焦失败
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
    
    % （3）一旦出现第二个Fdr的对比度最大，则步长减半
    while true
        
        IterateNumCurrent = IterateNumCurrent + 1;

        if  IterateNumCurrent>MaxIterateNum
            % 如果迭代次数超过最大迭代次数，则自聚焦失败
            FdrOfMaxContrast(m) = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent ;
            IsAutofocusSuccess(m) = false ;
            break ;
        end
        
        if FdrStepCurrent<=FdrStepFinal
            % 如果当前步长小于最终步长，则自聚焦成功，记录当前Fdr
            FdrOfMaxContrast(m) = FdrCurrent + (IndexOfMaxContrast-2)*FdrStepCurrent ;
            if FirstAutofocusSuccess==-1
                % 如果这是第一个自聚焦成功的聚焦深度，则记录下来
                FirstAutofocusSuccess = m ;
            end
            break ;
        end
        
        if IndexOfMaxContrast~=2
            % 如果第一个或第三个Fdr的对比度最大，则更新当前Fdr
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

%% 3、对自聚焦失败的聚焦深度进行处理
if(FirstAutofocusSuccess == -1)
    % 如果全部自聚焦失败，则使用初始值
    FdrOfMaxContrast = FdrInit ;
else
    % 如果没有全部自聚焦失败
    % 对于第一个聚焦深度，使用第一个自聚焦成功的Fdr值
    FdrOfMaxContrast(1) = FdrOfMaxContrast(FirstAutofocusSuccess);
    % 对于第二个到倒数第二个聚焦深度中自聚焦失败的聚焦深度，如果后一聚焦深度
    % 自聚焦成功则取前一个与后一个的均值，否则取前一个值。
    for k = 2 : NrNew/FdrEstimDepth-1
        if ~IsAutofocusSuccess(k)
            if IsAutofocusSuccess(k+1)
                FdrOfMaxContrast(k) = (FdrOfMaxContrast(k-1)+FdrOfMaxContrast(k+1))/2;
            else
                FdrOfMaxContrast(k) = FdrOfMaxContrast(k-1);
            end
        end
    end
    % 对于最后一个聚焦深度，如果自聚焦失败，则取前一个值
    k = k+1 ;
    if ~IsAutofocusSuccess(k)
        FdrOfMaxContrast(k) = FdrOfMaxContrast(k-1);
    end    
end

