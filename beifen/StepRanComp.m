%% 距离FFT，距离压缩，距离IFFT，去除距离向非完全积累点

Nwin = round(Nr*Br/Fs/2)*2 ;            %脉冲宽度内点数
WinRan = fftshift( [zeros((Nr-Nwin)/2,1) ; hamming(Nwin) ; zeros((Nr-Nwin)/2,1) ] ) ;       %窗函数
Hran = WinRan .* conj( fft( rectpuls((-Tp/2:1/Fs:Tp/2-1/Fs)',Tp) ...
    .*exp(SignOfChirpSlope*1j*pi*Br/Tp*(-Tp/2:1/Fs:Tp/2-1/Fs)'.^2) , Nr ) );                %匹配信号

FidReadReal = fopen( [FolderImageOutPut 'EchoReal.dat'] , 'r' ) ;
FidReadImag = fopen( [FolderImageOutPut 'EchoImag.dat'] , 'r' ) ;
FidWriteReal = fopen( 'RanCompReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'RanCompImag.dat' , 'w' ) ;

for m = 1 : Na
    
    Data = fread( FidReadReal , Nr , 'float32' ) + 1j * fread( FidReadImag , Nr , 'float32' ) ;
    Data = real(Data) - mean(real(Data)) + 1j*(imag(Data) - mean(imag(Data)));      %去直流
    Data = ifft( fft( Data ) .* Hran ) ;     %距离脉冲压缩+距离走动校正
    fwrite( FidWriteReal ,  real(Data(NrNewStart:(NrNewStart+NrNew-1))), 'float32' ) ;
    fwrite( FidWriteImag , imag(Data(NrNewStart:(NrNewStart+NrNew-1))) , 'float32' ) ;

end

fclose all ;
