%% 去除方位向高频分量

FidReadReal = fopen( 'RanCompTurnReal.dat' , 'r' ) ;
FidReadImag = fopen( 'RanCompTurnImag.dat' , 'r' ) ;
N_a = 3000;
% N_a = 2*round(2*(2*delta_theta*V/Lambda)/(Prf/Na)/2);
FidWriteReal = fopen( 'RanCompTurnReal1.dat' , 'w' ) ;
FidWriteImag = fopen( 'RanCompTurnImag1.dat' , 'w' ) ;
win = fftshift([zeros((Na - N_a)/2,1);hamming(N_a);zeros((Na - N_a)/2,1)]);
for m = 1 : NrNew
    
    Data = fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' ) ;
    Data = ifft( fft( Data ).*win ) ;
    fwrite( FidWriteReal , real(Data) , 'float32' ) ;
    fwrite( FidWriteImag , imag(Data) , 'float32' ) ;
    
end

fclose all ;