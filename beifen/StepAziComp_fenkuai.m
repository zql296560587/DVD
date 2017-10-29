%% 方位去斜，方位FFT

Ba = 2*delta_theta*V*cosd(ThetaVRefEstim)/Lambda;
% % N = round(Na*Prf/Ba/2)*2;
% N = 32768;
N = 51200;          %L波段FFT点数
N_Ba = round(Ba*N/Prf/2)*2-3600;
% 
% AlphaCalcu = 90;
% PixelAziNum = Na;
% GamaV = 0;
% PixelAziSpace = V/Prf;
% NumOfSubAzi = round(Na/4096);
% 
% DeltaAlpha = PixelAziNum*PixelAziSpace/sqrt(Rref^2-H^2)/pi*180;
% AlphaCalcuSubAzi = AlphaCalcu + (-1/2+(((NumOfSubAzi-1):-1:0)+1/2)/NumOfSubAzi)*DeltaAlpha;
% 
% trNew = 2*RnNew/C;
% TempSinBeta = H./(trNew*C/2);
% TempCosBeta = sqrt(1-TempSinBeta.^2);
% FdcSubAzi = 2*V/Lambda*cosd(GamaV)*TempCosBeta(NrNew/2)*(cosd(AlphaCalcuSubAzi)-cosd(AlphaCalcu));
% BoundsOfSubAzi = round((FdcSubAzi(2:NumOfSubAzi)+FdcSubAzi(1:NumOfSubAzi-1))/2/(Prf/Na))+Na/2+1;
% nStartSubAzi = [1 BoundsOfSubAzi];
% nEndSubAzi = [BoundsOfSubAzi+1 Na];
% 
% FdrEstimSubAzi = zeros(NrNew,NumOfSubAzi);
% TempA = -2*V^2/Lambda./(trNew*C/2)*sind(GamaV).*TempSinBeta*cosd(GamaV).*TempCosBeta;
% TempB = -2*V^2/Lambda./(trNew*C/2)*cosd(GamaV)^2.*TempCosBeta.^2;
% TempC = 0;
% for k = 1:NumOfSubAzi
%     DeltaFdrSubAzi = TempA*(cosd(AlphaCalcu)-cosd(AlphaCalcuSubAzi(k))) ...
%         + TempB * (cosd(AlphaCalcu)^2-cosd(AlphaCalcuSubAzi(k))^2);
%     FdrEstimSubAzi(:,k) = FdrEstim + DeltaFdrSubAzi;
% end

FidReadReal = fopen( 'AziIFftReal.dat' , 'r' ) ;
FidReadImag = fopen( 'AziIFftImag.dat' , 'r' ) ;

FidWriteReal = fopen( 'AziCompReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'AziCompImag.dat' , 'w' ) ;

WinAzi = hamming(Na) ;

% win = fftshift([zeros(1000,1);hamming(Na-2000);zeros(1000,1)]);

for m = 1 : NrNew
    Data =  fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' );
%     Data = ifft(fft(Data).*win);
%% specan 
%     dataComp = zeros(Na,1);
%     for n = 1:NumOfSubAzi
        Hazi = WinAzi.*exp( -2j*pi*(FdcEstim(m)-FdcRefEstim)*ta -1j*pi*FdrEstim(m)*ta.^2 ...
            - 1j*pi*( -2*V^3*cosd(ThetaVEstim(m))*sind(ThetaVEstim(m))^2/(3*Lambda*RnNew(m)^2)) *ta.^3 );
%             - 1j*pi*( V^4*sind(ThetaVEstim(m))^2*(5*cosd(ThetaVEstim(m))^2-1)/8/Lambda/RnNew(m)^3)*ta.^4 );
        dataSub = Data.*Hazi;
        dataComp = fftshift(fft(fftshift(dataSub),N));
%         dataComp(nStartSubAzi(n):nEndSubAzi(n)) = dataSub(nStartSubAzi(n)+N/2-Na/2:nEndSubAzi(n)+N/2-Na/2);
%     end

    fwrite( FidWriteReal , real(dataComp(Na/2-N_Ba/2 : N_Ba/2+Na/2-1)) , 'float32' ) ;
    fwrite( FidWriteImag , imag(dataComp(Na/2-N_Ba/2 : N_Ba/2+Na/2-1)) , 'float32' ) ;
    
end
% for m = 1 : NrNew
%     Data =  fread( FidReadReal , Na , 'float32' ) + 1j * fread( FidReadImag , Na , 'float32' );
%     Data = ifft(fft(Data).*win);
% %% specan 
%     dataComp = zeros(Na,1);
%     for n = 1:NumOfSubAzi
%         Hazi = WinAzi.*exp( -2j*pi*(FdcEstim(m)-FdcRefEstim)*ta -1j*pi*FdrEstimSubAzi(m,n)*ta.^2 ...
%             - 1j*pi*( -2*V^3*cosd(ThetaVEstim(m))*sind(ThetaVEstim(m))^2/(3*Lambda*RnNew(m)^2)) *ta.^3 );
% %             - 1j*pi*( V^4*sind(ThetaVEstim(m))^2*(5*cosd(ThetaVEstim(m))^2-1)/8/Lambda/RnNew(m)^3)*ta.^4 );
%         dataSub = Data.*Hazi;
%         dataSub = fftshift(fft(fftshift(dataSub),N));
%         dataComp(nStartSubAzi(n):nEndSubAzi(n)) = dataSub(nStartSubAzi(n)+N/2-Na/2:nEndSubAzi(n)+N/2-Na/2);
%     end
% 
%     fwrite( FidWriteReal , real(dataComp(Na/2-N_Ba/2 : N_Ba/2+Na/2-1)) , 'float32' ) ;
%     fwrite( FidWriteImag , imag(dataComp(Na/2-N_Ba/2 : N_Ba/2+Na/2-1)) , 'float32' ) ;
%     
% end

fclose all ;

