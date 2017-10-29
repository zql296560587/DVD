%% ����FFT������ѹ��������IFFT��ȥ�����������ȫ���۵�

Nwin = round(Nr*Br/Fs/2)*2 ;            %�������ڵ���
WinRan = fftshift( [zeros((Nr-Nwin)/2,1) ; hamming(Nwin) ; zeros((Nr-Nwin)/2,1) ] ) ;       %������
Hran = WinRan .* conj( fft( rectpuls((-Tp/2:1/Fs:Tp/2-1/Fs)',Tp) ...
    .*exp(SignOfChirpSlope*1j*pi*Br/Tp*(-Tp/2:1/Fs:Tp/2-1/Fs)'.^2) , Nr ) );                %ƥ���ź�

FidReadReal = fopen( [FolderImageOutPut 'EchoReal.dat'] , 'r' ) ;
FidReadImag = fopen( [FolderImageOutPut 'EchoImag.dat'] , 'r' ) ;
FidWriteReal = fopen( 'RanCompReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'RanCompImag.dat' , 'w' ) ;

for m = 1 : Na
    
    Data = fread( FidReadReal , Nr , 'float32' ) + 1j * fread( FidReadImag , Nr , 'float32' ) ;
    Data = real(Data) - mean(real(Data)) + 1j*(imag(Data) - mean(imag(Data)));      %ȥֱ��
    Data = ifft( fft( Data ) .* Hran ) ;     %��������ѹ��+�����߶�У��
    fwrite( FidWriteReal ,  real(Data(NrNewStart:(NrNewStart+NrNew-1))), 'float32' ) ;
    fwrite( FidWriteImag , imag(Data(NrNewStart:(NrNewStart+NrNew-1))) , 'float32' ) ;

end

fclose all ;
