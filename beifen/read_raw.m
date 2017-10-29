function read_raw(inputfile,outputfilereal,outputfileimag,Na,Nr,index,loop)
%从data文件中读取当前孔径的回波数据
Nr = 2*Nr;
fd1 = fopen(inputfile,'rb');
fd2 = fopen(outputfilereal,'wb');
fd3 = fopen(outputfileimag,'wb');
data1 = zeros(2*Nr,1);
fseek(fd1,0,'eof');
N = ftell(fd1);
num = N/(Nr*4*2+128*8)/16384;
fseek(fd1,0,'bof');
% fseek(fd1,16384*2*(Nr*4*2+128*8),'bof');
fseek(fd1,32768*(Nr*4+128*8)*35+(Nr*4+128*8)*8192*(loop-1),'bof');     %L 
% fseek(fd1,32768*(Nr*4+128*8)*35,'bof');
Ass1 = [];
Ass2 = [];
Ass3 = [];
Ass4 = [];
% a = 16391;
a = 32783;
% echo = [];
if index == 1
    for i = 1:Na
        dataI = [];
        dataQ = [];
        fseek(fd1,256*2,'cof');
%         Ass = [Ass,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
        dataI = [dataI;fread(fd1,128,'int8')];
%         Ass = [Ass,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
        dataQ = [dataQ;fread(fd1,128,'int8')];
        for n = 1:(Nr-256)/256
            fseek(fd1,256*2,'cof');
            dataI = [dataI;fread(fd1,256,'int8')];
            dataQ = [dataQ;fread(fd1,256,'int8')];
        end
        fseek(fd1,256*2,'cof');
%         fseek(fd1,128,'cof');
        dataI = [dataI;fread(fd1,128,'int8')];     %实部
%         Ass1 = [Ass1,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
%         fseek(fd1,128,'cof');
        dataQ = [dataQ;fread(fd1,128,'int8')];     %虚部
        fseek(fd1,128,'cof');
%         Ass1 = [Ass1,fread(fd1,128,'uint8')];
        data1(1:2:end) = dataI;
        data1(2:2:end) = dataQ;
    %     if(i == 16129)
    %    data(1:2:end) = dataI;
    %     end
        i
    %     echo = [echo,dataI+1j*dataQ];
        fwrite(fd2,dataI(1:Nr/2),'float32');
        fwrite(fd3,dataQ(1:Nr/2),'float32');
    end
elseif index == 5
    for i = 1:Na
        if i==1200
            1;
        end
        dataI = [];
        dataQ = [];
        Ass1 = [Ass1,fread(fd1,128,'uint8')];
        num_prf = Ass1(8,i)*256^3+Ass1(7,i)*256^2+Ass1(6,i)*256+Ass1(5,i);
        a = num_prf-a;
%         fseek(fd1,128,'cof');
        dataI = [dataI;fread(fd1,64,'int16')];
        Ass2 = [Ass2,fread(fd1,128,'uint8')];
%         fseek(fd1,128,'cof');
        dataQ = [dataQ;fread(fd1,64,'int16')]; 
        Ass3 = [Ass3,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
        Ass4 = [Ass4,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
%         Ass = [Ass,fread(fd1,128,'uint8')];
        for n = 1:(Nr-128)/128
            dataI = [dataI;fread(fd1,128,'int16')];
            dataQ = [dataQ;fread(fd1,128,'int16')];
            fseek(fd1,256*2,'cof');
        end
        dataI = [dataI;fread(fd1,64,'int16')];     %实部
%         Ass1 = [Ass1,fread(fd1,128,'uint8')];
        fseek(fd1,128,'cof');
        dataQ = [dataQ;fread(fd1,64,'int16')];     %虚部
        fseek(fd1,128,'cof');
        fseek(fd1,256*2,'cof');
        i
        if a == 1
            fwrite(fd2,dataI(1:15000),'float32');
            fwrite(fd3,dataQ(1:15000),'float32');
            dataI1 = dataI;
            dataQ1 = dataQ;
        else
            fwrite(fd2,(dataI1(1:15000)+dataI(1:15000))/2,'float32');
            fwrite(fd3,(dataQ1(1:15000)+dataQ(1:15000))/2,'float32');
            fwrite(fd2,dataI(1:15000),'float32');
            fwrite(fd3,dataQ(1:15000),'float32');           
        end
        a = num_prf;
    end
end
        
Nr = Nr/2;        %L
fclose all;

