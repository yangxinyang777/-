   %空间后方交会迭代法
   %像点坐标以mm为单位
    %*************************************************************
    clc
    clear
    format long
    fid=fopen('相机参数.txt','r');
    A=fscanf(fid,'%f');
    fp=fopen('外方位元素初值.txt','r');
    B=fscanf(fp,'%f');
    fm=fopen('控制点文件.txt','r');
    C=fscanf(fm,'%f');
    fr=fopen('.txt','w');%修改为自己的文件路径和名称
    x0=A(1);
    y0=A(2);%像主点坐标
    f=A(3);%主距
    m=A(4);%CCD像素大小
    k1=A(5);
    k2=A(6);%相机镜头的径向畸变差参数
    p1=A(7);
    p2=A(8);%切向畸变差参数
    a=A(9);
    b=A(10);%像元仿射变换参数
    Xsl=B(1);
    Ysl=B(2);
    Zsl=B(3);
    fil=B(4);
    wal=B(5);
    kal=B(6);%影像外方位元素
    n=C(1);%控制点个数
    for i=1:n
        PTNUM(i)=C(6*(i-1)+2);
        X(i)=C(6*(i-1)+3);
        Y(i)=C(6*(i-1)+4);
        Z(i)=C(6*(i-1)+5);%控制点坐标，右手坐标系
        xl(i)=C(6*(i-1)+6)*m;
        yl(i)=C(6*(i-1)+7)*m;%左影像像点坐标
    end
    %单像后方交会
    fprintf(fr,'单像后方交会结果输出\r\n左外方位元素改正值\r\n           dXl           dYl          dZl          dfil        dwal         dkal\r\n');
    dfil=0.1;%外方位元素角元素初值
    dwal=0.1;
    dkal=0.1;
    sigmal=0.1;
    count0=0;
    while abs(dfil)>1e-7||abs(dwal)>1e-7||abs(dkal)>1e-7%||sigmal>m
        count0=count0+1;
        %计算旋转矩阵R 
        al1=cos(fil)*cos(kal)-sin(fil)*sin(wal)*sin(kal);
        al2=-cos(fil)*sin(kal)-sin(fil)*sin(wal)*cos(kal);
        al3=-sin(fil)*cos(wal);
        bl1=cos(wal)*sin(kal);
        bl2=cos(wal)*cos(kal);
        bl3=-sin(wal);
        cl1=sin(fil)*cos(kal)+cos(fil)*sin(wal)*sin(kal);
        cl2=-sin(fil)*sin(kal)+cos(fil)*sin(wal)*cos(kal);
        cl3=cos(fil)*cos(wal);                               %旋转矩阵R
    for i=1:n
        rl=sqrt((xl(i)-x0)*(xl(i)-x0)+(yl(i)-y0)*(yl(i)-y0));
        dxl(i)=(xl(i)-x0)*(k1*rl*rl+k2*rl*rl*rl*rl)+p1*(rl*rl+2*(xl(i)-x0)*(xl(i)-x0))+2*p2*(xl(i)-x0)*(yl(i)-y0)+a*(xl(i)-x0)+b*(yl(i)-y0);
        dyl(i)=(yl(i)-y0)*(k1*rl*rl+k2*rl*rl*rl*rl)+p2*(rl*rl+2*(yl(i)-y0)*(yl(i)-y0))+2*p1*(xl(i)-x0)*(yl(i)-y0);                           %畸变差
        Xhl(i)=al1*(X(i)-Xsl)+bl1*(Y(i)-Ysl)+cl1*(Z(i)-Zsl);
        Yhl(i)=al2*(X(i)-Xsl)+bl2*(Y(i)-Ysl)+cl2*(Z(i)-Zsl);
        Zhl(i)=al3*(X(i)-Xsl)+bl3*(Y(i)-Ysl)+cl3*(Z(i)-Zsl);%共线方程分子、分母
        xlj(i)=-f*Xhl(i)/Zhl(i);
        ylj(i)=-f*Yhl(i)/Zhl(i);%控制点像点坐标的近似值
        %注意：共线方程处的x=像平面坐标x-x0（内方位元素）-相片畸变差改正dx
%       dxl(i)=0;
%       dyl(i)=0;
        ddxl(i)=xl(i)-x0-dxl(i);
        ddyl(i)=yl(i)-y0-dyl(i);
        al11(i)=(al1*f+al3*ddxl(i))/Zhl(i);
        al12(i)=(bl1*f+bl3*ddxl(i))/Zhl(i);
        al13(i)=(cl1*f+cl3*ddxl(i))/Zhl(i);
        al14(i)=ddyl(i)*sin(wal)-(ddxl(i)*(ddxl(i)*cos(kal)-ddyl(i)*sin(kal))/f+f*cos(kal))*cos(wal);
        al15(i)=-f*sin(kal)-ddxl(i)*(ddxl(i)*sin(kal)+ddyl(i)*cos(kal))/f;
        al16(i)=ddyl(i);
        al21(i)=(al2*f+al3*ddyl(i))/Zhl(i);
        al22(i)=(bl2*f+bl3*ddyl(i))/Zhl(i);
        al23(i)=(cl2*f+cl3*ddyl(i))/Zhl(i);
        al24(i)=-ddxl(i)*sin(wal)-(ddyl(i)*(ddxl(i)*cos(kal)-ddyl(i)*sin(kal))/f-f*sin(kal))*cos(wal);
        al25(i)=-f*cos(kal)-ddyl(i)*(ddxl(i)*sin(kal)+ddyl(i)*cos(kal))/f;
        al26(i)=-ddxl(i);
        Al(2*i-1,1)=al11(i);
        Al(2*i,1)=al21(i);
        Al(2*i-1,2)=al12(i);
        Al(2*i,2)=al22(i);
        Al(2*i-1,3)=al13(i);
        Al(2*i,3)=al23(i);
        Al(2*i-1,4)=al14(i);
        Al(2*i,4)=al24(i);
        Al(2*i-1,5)=al15(i);
        Al(2*i,5)=al25(i);
        Al(2*i-1,6)=al16(i);
        Al(2*i,6)=al26(i);%A矩阵
        lxl(i)=xl(i)-x0-dxl(i)-xlj(i);
        lyl(i)=yl(i)-y0-dyl(i)-ylj(i);
        Ll(2*i-1)=lxl(i);
        Ll(2*i)=lyl(i);%L矩阵    
    end
    dXl=inv(Al'*Al)*Al'*Ll';
    Xsl=Xsl+dXl(1);
    Ysl=Ysl+dXl(2);
    Zsl=Zsl+dXl(3);
    fil=fil+dXl(4);
    wal=wal+dXl(5);
    kal=kal+dXl(6);
    dfil=dXl(4);
    dwal=dXl(5);
    dkal=dXl(6);
    Vl=Al*dXl-Ll';
    sigmal=abs(sqrt(Vl'*Vl/(2*n-6)));
    fprintf(fr,'第%d次迭代 %f    %f    %f      %f    %f    %f\r\n',count0,dXl);
    end
    plot(dXl);
    title('后方交会改正数迭代数图');
    xlabel('迭代次数');
    ylabel('dXl');
    grid on

    fprintf(fr,'外方元素值\r\n  Xsl             Ysl            Zsl          fil           wal          kal\r\n');
    fprintf(fr,'%f       %f     %f    %f    %f     %f\r\n',Xsl,Ysl,Zsl,fil,wal,kal);
    fprintf(fr,'单位权中误差（单位:像素）%f\r\n\r\n',sigmal/m);
    fprintf(fr,'*********************************************分格线***************************************************');