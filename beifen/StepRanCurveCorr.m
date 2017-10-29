%% 距离FFT，距离弯曲校正，距离IFFT
%% 注释为分块校正，目前距离向点数不大，可不用分块
% Num_Rcmc = round(NrNew/2048);
% Rcmc = zeros(Na,Num_Rcmc);
% 
% BoundsOfRcmc = round(((1:Num_Rcmc-1)/Num_Rcmc)*NrNew);
% nStartSub = [1 BoundsOfRcmc];
% nEndSub = [BoundsOfRcmc+1 NrNew];
% NSub = round((nStartSub + nEndSub)./2);
% 
% for i = 1:Num_Rcmc
%     ThetaVRefEstim0 = asind(Lambda*FdcEstim(NSub(i))/2/V);
%     for m = 1:Na
%         if 1-(Lambda*(fa(m)+FdcEstim(NSub(i)))/2/V)^2>0
%              Rcmc(m,i) = (RminNew + bin_r*NSub(i))*cosd(ThetaVRefEstim0)*( (1-Lambda*(fa(m)+FdcEstim(NSub(i)))/2/V*sind(ThetaVRefEstim0)) ...
%                     ./sqrt(1-(Lambda*(fa(m)+FdcEstim(NSub(i)))/2/V).^2)-cosd(ThetaVRefEstim0) ) ; % 距离弯曲校正量
%         end
%     end
% end
Rcmc = zeros(Na,1);
for m = 1:Na
    if 1-(Lambda*(fa(m)+FdcRefEstim)/2/V)^2>0
         Rcmc(m) = Rref*cosd(ThetaVRefEstim) * ( (1-Lambda*(fa(m)+FdcRefEstim)/2/V*sind(ThetaVRefEstim)) ...
                ./sqrt(1-(Lambda*(fa(m)+FdcRefEstim)/2/V).^2)-cosd(ThetaVRefEstim) ) ; % 距离弯曲校正量
    end
end
% Rcmc = zeros(Na,1);
FidReadReal = fopen( 'AziFftTurnReal.dat' , 'r' ) ;
FidReadImag = fopen( 'AziFftTurnImag.dat' , 'r' ) ;

FidWriteReal = fopen( 'RanCurveCorrReal.dat' , 'w' ) ;
FidWriteImag = fopen( 'RanCurveCorrImag.dat' , 'w' ) ;

for m = 1 : Na

    Data = fread( FidReadReal , NrNew , 'float32' ) + 1j * fread( FidReadImag , NrNew , 'float32' ) ;
    Data = ifft( fft(Data) .* exp(2j*pi*2*Rcmc(m)/C*frNew )) ;
    
    fwrite( FidWriteReal , real(Data) , 'float32' ) ;
    fwrite( FidWriteImag , imag(Data) , 'float32' ) ;

end
% for m = 1 : Na
%     Data0 = zeros(NrNew,1);
%     Data = fread( FidReadReal , NrNew , 'float32' ) + 1j * fread( FidReadImag , NrNew , 'float32' ) ;
%     for i = 1:Num_Rcmc
%         Data1 = ifft( fft(Data) .* exp(2j*pi*2*Rcmc(m,i)/C*frNew ) ) ;
%         Data0( nStartSub(i) : nEndSub(i)) = Data1( nStartSub(i) : nEndSub(i) );
%     end
%     fwrite( FidWriteReal , real(Data) , 'float32' ) ;
%     fwrite( FidWriteImag , imag(Data) , 'float32' ) ;
% 
% end
%             .*exp(-1j*pi*2*Rref*cosd(ThetaVRefEstim)*Lambda/C^2*(Lambda.*fa(m)/2/V)^2 ...
%             /(sqrt(1-(Lambda*(fa(m)+FdcRefEstim)/2/V)^2))^3*fr.^2));
fclose all ;
