   %�ռ�󷽽��������
   %���������mmΪ��λ
    %*************************************************************
    clc
    clear
    format long
    fid=fopen('�������.txt','r');
    A=fscanf(fid,'%f');
    fp=fopen('�ⷽλԪ�س�ֵ.txt','r');
    B=fscanf(fp,'%f');
    fm=fopen('���Ƶ��ļ�.txt','r');
    C=fscanf(fm,'%f');
    fr=fopen('.txt','w');%�޸�Ϊ�Լ����ļ�·��������
    x0=A(1);
    y0=A(2);%����������
    f=A(3);%����
    m=A(4);%CCD���ش�С
    k1=A(5);
    k2=A(6);%�����ͷ�ľ����������
    p1=A(7);
    p2=A(8);%�����������
    a=A(9);
    b=A(10);%��Ԫ����任����
    Xsl=B(1);
    Ysl=B(2);
    Zsl=B(3);
    fil=B(4);
    wal=B(5);
    kal=B(6);%Ӱ���ⷽλԪ��
    n=C(1);%���Ƶ����
    for i=1:n
        PTNUM(i)=C(6*(i-1)+2);
        X(i)=C(6*(i-1)+3);
        Y(i)=C(6*(i-1)+4);
        Z(i)=C(6*(i-1)+5);%���Ƶ����꣬��������ϵ
        xl(i)=C(6*(i-1)+6)*m;
        yl(i)=C(6*(i-1)+7)*m;%��Ӱ���������
    end
    %����󷽽���
    fprintf(fr,'����󷽽��������\r\n���ⷽλԪ�ظ���ֵ\r\n           dXl           dYl          dZl          dfil        dwal         dkal\r\n');
    dfil=0.1;%�ⷽλԪ�ؽ�Ԫ�س�ֵ
    dwal=0.1;
    dkal=0.1;
    sigmal=0.1;
    count0=0;
    while abs(dfil)>1e-7||abs(dwal)>1e-7||abs(dkal)>1e-7%||sigmal>m
        count0=count0+1;
        %������ת����R 
        al1=cos(fil)*cos(kal)-sin(fil)*sin(wal)*sin(kal);
        al2=-cos(fil)*sin(kal)-sin(fil)*sin(wal)*cos(kal);
        al3=-sin(fil)*cos(wal);
        bl1=cos(wal)*sin(kal);
        bl2=cos(wal)*cos(kal);
        bl3=-sin(wal);
        cl1=sin(fil)*cos(kal)+cos(fil)*sin(wal)*sin(kal);
        cl2=-sin(fil)*sin(kal)+cos(fil)*sin(wal)*cos(kal);
        cl3=cos(fil)*cos(wal);                               %��ת����R
    for i=1:n
        rl=sqrt((xl(i)-x0)*(xl(i)-x0)+(yl(i)-y0)*(yl(i)-y0));
        dxl(i)=(xl(i)-x0)*(k1*rl*rl+k2*rl*rl*rl*rl)+p1*(rl*rl+2*(xl(i)-x0)*(xl(i)-x0))+2*p2*(xl(i)-x0)*(yl(i)-y0)+a*(xl(i)-x0)+b*(yl(i)-y0);
        dyl(i)=(yl(i)-y0)*(k1*rl*rl+k2*rl*rl*rl*rl)+p2*(rl*rl+2*(yl(i)-y0)*(yl(i)-y0))+2*p1*(xl(i)-x0)*(yl(i)-y0);                           %�����
        Xhl(i)=al1*(X(i)-Xsl)+bl1*(Y(i)-Ysl)+cl1*(Z(i)-Zsl);
        Yhl(i)=al2*(X(i)-Xsl)+bl2*(Y(i)-Ysl)+cl2*(Z(i)-Zsl);
        Zhl(i)=al3*(X(i)-Xsl)+bl3*(Y(i)-Ysl)+cl3*(Z(i)-Zsl);%���߷��̷��ӡ���ĸ
        xlj(i)=-f*Xhl(i)/Zhl(i);
        ylj(i)=-f*Yhl(i)/Zhl(i);%���Ƶ��������Ľ���ֵ
        %ע�⣺���߷��̴���x=��ƽ������x-x0���ڷ�λԪ�أ�-��Ƭ��������dx
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
        Al(2*i,6)=al26(i);%A����
        lxl(i)=xl(i)-x0-dxl(i)-xlj(i);
        lyl(i)=yl(i)-y0-dyl(i)-ylj(i);
        Ll(2*i-1)=lxl(i);
        Ll(2*i)=lyl(i);%L����    
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
    fprintf(fr,'��%d�ε��� %f    %f    %f      %f    %f    %f\r\n',count0,dXl);
    end
    plot(dXl);
    title('�󷽽��������������ͼ');
    xlabel('��������');
    ylabel('dXl');
    grid on

    fprintf(fr,'�ⷽԪ��ֵ\r\n  Xsl             Ysl            Zsl          fil           wal          kal\r\n');
    fprintf(fr,'%f       %f     %f    %f    %f     %f\r\n',Xsl,Ysl,Zsl,fil,wal,kal);
    fprintf(fr,'��λȨ������λ:���أ�%f\r\n\r\n',sigmal/m);
    fprintf(fr,'*********************************************�ָ���***************************************************');