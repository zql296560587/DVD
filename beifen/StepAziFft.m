%% ·½Î»FFT

FidReadReal = fopen( 'ResiRanWalkCorrTurnReal.dat' , 'r' ) ;
FidReadImag = fopen( 'ResiRanWalkCorrTurnImag.dat' , 'r' ) ;

FidWriteReal = fopen( 'AziFftReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'AziFftImag.dat' , 'w' ) ;

for m = 1 : NrNew
    
    Data = fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' ) ;
    Data = ( fft( Data ) ) ;
    fwrite( FidWriteReal , real(Data) , 'float32' ) ;
    fwrite( FidWriteImag , imag(Data) , 'float32' ) ;
    
end

fclose all ;
