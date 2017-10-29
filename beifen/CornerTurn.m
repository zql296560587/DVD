function CornerTurn( FileIn , FileOut , RowNum , ColNum , DivNum )

%Ϊ�������ݾ���ת��ʱ���������ռ���ڴ棬����Ҫ�����ݾ���ֳɼ��ݷֱ�ת

FidIn = fopen( FileIn , 'r' ) ;
FidOut = fopen( FileOut , 'w' ) ;

%���з����ݾ���
SubColNum = ColNum / DivNum ; %��תʱÿһ�ݵ��еĴ�С

Data = zeros( RowNum , SubColNum ) ; %��ʱ���ÿһС�ݵ����ݾ���

for k = 1:DivNum
    
    frewind( FidIn );
    SpaceByteStart = 4*(k-1)*SubColNum ;
    SpaceByteEnd = 4*(DivNum-k)*SubColNum;
    
    for m = 1 : RowNum
        fseek( FidIn , SpaceByteStart , 'cof' ) ;
        Temp = fread( FidIn , SubColNum , 'float32' ) ; %��ָ��ָ��ÿһ�е�Ҫ������е���ʼ��
        Data(m,:) = Temp' ;
        fseek( FidIn , SpaceByteEnd , 'cof' ) ; %��ָ��ָ��ÿһ�еĿ�ʼ
    end
    
    for n = 1 : SubColNum
        fwrite( FidOut , Data(:,n) , 'float32' ) ; %�Ӿ���ķ�ת
    end
    
end

fclose(FidIn);
fclose(FidOut);

