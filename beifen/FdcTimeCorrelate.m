%% 时域相关法估计多普勒中心

FdcEstimDepth = 1 ;    % 每FdcEstimDepth个距离门估计一次，要求能被距离门个数NrNew整除

FidReadReal = fopen( 'RanCompTurnReal1.dat' , 'r' ) ;
FidReadImag = fopen( 'RanCompTurnImag1.dat' , 'r' ) ;

FdcOfTimeCorr = zeros( NrNew/FdcEstimDepth , 1 ) ;
% fseek(FidReadReal,(polyfit_fdc_start-1)*Na*4,'bof');
% fseek(FidReadImag,(polyfit_fdc_start-1)*Na*4,'bof');
for m = 1 : NrNew/FdcEstimDepth
    
    Data1 = fread( FidReadReal , Na , 'float32' ) + 1i * fread( FidReadImag , Na , 'float32' ) ;
    
    FdcOfTimeCorr(m) = angle( sum( Data1(2:end).*conj(Data1(1:end-1)) ) ) / (2*pi) * Prf ;          %计算fdc值  
    
    fseek( FidReadReal , (FdcEstimDepth-1)*Na*4 , 'cof' ) ;
    fseek( FidReadImag , (FdcEstimDepth-1)*Na*4 , 'cof' ) ;
    
end

fclose all ;
