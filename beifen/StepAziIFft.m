%% ∑ΩŒªœÚIFFT

FidReadReal = fopen( 'RanCurveCorrTurnReal.dat' , 'r' ) ;
FidReadImag = fopen( 'RanCurveCorrTurnImag.dat' , 'r' ) ;
FidWriteReal = fopen( 'AziIFftReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'AziIFftImag.dat' , 'w' ) ;

for m = 1 : NrNew
    
    Data = fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' ) ;
    Data = ifft( ( Data ) ) ;
    fwrite( FidWriteReal , real(Data) , 'float32' ) ;
    fwrite( FidWriteImag , imag(Data) , 'float32' ) ;
    
end

fclose all ;