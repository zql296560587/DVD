%% 距离FFT，残余距离走动校正，去残余多普勒中心

FidReadReal = fopen( 'RanCompReal1.dat' , 'r' ) ;
FidReadImag = fopen( 'RanCompImag1.dat' , 'r' ) ;
FidWriteReal = fopen( 'ResiRanWalkCorrReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'ResiRanWalkCorrImag.dat' , 'w' ) ;

for m = 1 : Na
    
    Data = fread( FidReadReal , NrNew , 'float32' ) + 1j * fread( FidReadImag , NrNew , 'float32' ) ;
    
    Data = fft(Data) .* exp( -2j*pi*(FdcRefEstim-FdcRefCalcu)*ta(m)*(frNew*Lambda/C)-2j*pi*FdcRefEstim*ta(m) ) ;     
               %残余距离走动校正，去除多普勒中心和弯曲相位补偿
    
    Data = ifft(Data) ;
    
    fwrite( FidWriteReal , real(Data) , 'float32' ) ;
    fwrite( FidWriteImag , imag(Data) , 'float32' ) ;
    
end

fclose all ;

   
    