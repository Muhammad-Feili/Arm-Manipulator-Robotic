function proj
    function H=A(th, d, a, alpha)
             H=[cos(th),   -sin(th)*cos(alpha),  sin(th)*sin(alpha),   a*cos(th);  
                sin(th),   cos(th)*cos(alpha),   -cos(th)*sin(alpha),  a*sin(th);
                0,         sin(alpha),           cos(alpha),           d        ;
                0,         0,                    0,                    1       ];
    end

    function H=FK(th1,th2,th3,th4)
        C1=cos(th1); C2=cos(th2); C3=cos(th3); C34=cos(th3+th4);
        S1=sin(th1); S2=sin(th2); S3=sin(th3); S34=sin(th3+th4);
        
        L1=6.8; L2=9.3; L3=10.9; L4=8.5;
        
        nx=C1*C2*C34+S1*S34;
        ny=S1*C2*C34-C1*S34;
        nz=S2*C34;
        ox=C1*S2;
        oy=S1*S2;
        oz=-C2;
        ax=C1*C2*S34-S1*C34;
        ay=S1*C2*S34+C1*C34;
        az=S2*S34;
        px=L4*C1*C2*C34+L3*C1*C2*C3+L4*S1*S34+L3*S1*S3+L2*C1*C2;
        py=L4*S1*C2*C34+L3*S1*C2*C3-L4*C1*S34-L3*C1*S3+L2*S1*C2;
        pz=L4*S2*C34+L3*S2*C3+L2*S2+L1;
        
        H=[ nx ox ax px;
            ny oy ay py;
            nz oz az pz;
            0  0  0  1];
    end
    function H=IK(X,Y,Z)
        
        L1=6.8; L2=9.3; L3=10.9;
        
        th1=atan(Y/X);
        %S1=sin(th1); C1=cos(th1);
        
        TH1=th1+pi;
        
        R=sqrt(X^2+Y^2);
        C3=(R^2+(Z-L1)^2-L2^2-L3^2)/(2*L2*L3);
        S3=sqrt(1-C3^2);
        th3=atan(S3/C3);
        
        TH3=th3+pi;
        
        C2=L2+L3*C3;
        S2=L3*S3;
        th2=atan(S2/C2);
        
        TH2=atan((L3*sin(TH3))/(L2+L3*cos(TH3)));
        
        %C123=cos(th1+th2+th3); S123=sin(th1+th2+th3);
        %C23=cos(th2+th3); S23=sin(th2+th3);
        
        A11=A(th1,L1,0,pi/2);
        A12=A(TH2,0,L2,0);
        A13=A(th3,0,L2,0);
        
        T13=A11*A12*A13;
        R13=T13(1,3); R23=T13(2,3); R33=T13(3,3);

        s41=-1*(cos(th1+TH2+th3)*R13+cos(TH2+th3)*sin(th1)*R23+sin(TH2+th3)*R33);
        c41=-1*(sin(TH2+th3)*cos(th1)*R13+sin(th1+TH2+th3)*R23-cos(TH2+th3)*R33);
        th41=atan(s41/c41);
        
        A21=A(TH1,L1,0,pi/2);
        A22=A(TH2,0,L2,0);
        A23=A(th3,0,L2,0);
        
        T23=A21*A22*A23;
        R13=T23(1,3); R23=T23(2,3); R33=T23(3,3);
        
        s42=-1*(cos(TH1+TH2+th3)*R13+cos(TH2+th3)*sin(TH1)*R23+sin(TH2+th3)*R33);
        c42=-1*(sin(TH2+th3)*cos(TH1)*R13+sin(TH1+TH2+th3)*R23-cos(TH2+th3)*R33);
        th42=atan(s42/c42);
        
        A31=A(th1,L1,0,pi/2);
        A32=A(th2,0,L2,0);
        A33=A(TH3,0,L2,0);
        
        T33=A31*A32*A33;
        R13=T33(1,3); R23=T33(2,3); R33=T33(3,3);
        
        s43=-1*(cos(th1+th2+TH3)*R13+cos(th2+TH3)*sin(th1)*R23+sin(th2+TH3)*R33);
        c43=-1*(sin(th2+TH3)*cos(th1)*R13+sin(th1+th2+TH3)*R23-cos(th2+TH3)*R33);
        th43=atan(s43/c43);
        
        A41=A(TH1,L1,0,pi/2);
        A42=A(th2,0,L2,0);
        A43=A(TH3,0,L2,0);
        
        T43=A41*A42*A43;
        R13=T43(1,3); R23=T43(2,3); R33=T43(3,3);
        
        s44=-1*(cos(TH1+th2+TH3)*R13+cos(th2+TH3)*sin(TH1)*R23+sin(th2+TH3)*R33);
        c44=-1*(sin(th2+TH3)*cos(TH1)*R13+sin(TH1+th2+TH3)*R23-cos(th2+TH3)*R33);
        th44=atan(s44/c44);
        
        H=[th1 TH2 th3 th41;
           TH1 TH2 th3 th42;
           th1 th2 TH3 th43+pi;
           TH1 th2 TH3 th44+pi];
    end


f=figure(1);
set(gcf,'Position',[100 100 680 370]);
subplot(2,2,2);
%inverse kinematic button
inverse=uicontrol(f);
inverse.String='Inverse';
inverse.Units='pixels';
inverse.FontSize=12;
inverse.Position=[185 100 65 25];
inverse.Callback=@Inverse;
%forward kinematic button
inverse=uicontrol(f);
inverse.String='Forward';
inverse.Units='pixels';
inverse.FontSize=12;
inverse.Position=[65 100 65 25];
inverse.Callback=@Forward;
%X_position
PX=uicontrol(f);
PX.Style='edit';
PX.Position=[200 310 50 25];
PX.Units='pixels';
%X_label
LX=uicontrol(f);
LX.Style='text';
LX.Position=[160 310 40 25];
LX.String='X';
LX.FontSize=16;
%Y_position
PY=uicontrol(f);
PY.Style='edit';
PY.Position=[200 235 50 25];
PY.Units='pixels';
%Y_label
LY=uicontrol(f);
LY.Style='text';
LY.Position=[160 235 40 25];
LY.String='Y';
LY.FontSize=16;
%Z_position
PZ=uicontrol(f);
PZ.Style='edit';
PZ.Position=[200 160 50 25];
PZ.Units='pixels';
%Z_label
LZ=uicontrol(f);
LZ.Style='text';
LZ.Position=[160 160 40 25];
LZ.String='Z';
LZ.FontSize=16;

%joints variables

%theta1 position
Th1=uicontrol(f);
Th1.Style='edit';
Th1.Position=[80 310 50 25];
Th1.Units='pixels';
%theta1 label
LTH1=uicontrol(f);
LTH1.Style='text';
LTH1.Position=[30 310 40 25];
LTH1.String='TH1';
LTH1.FontSize=16;
%theta2 position
Th2=uicontrol(f);
Th2.Style='edit';
Th2.Position=[80 260 50 25];
Th2.Units='pixels';
%theta2 label
LTH2=uicontrol(f);
LTH2.Style='text';
LTH2.Position=[30 260 40 25];
LTH2.String='TH2';
LTH2.FontSize=16;
%theta3 position
Th3=uicontrol(f);
Th3.Style='edit';
Th3.Position=[80 210 50 25];
Th3.Units='pixels';
%theta3 label
LTH3=uicontrol(f);
LTH3.Style='text';
LTH3.Position=[30 210 40 25];
LTH3.String='TH3';
LTH3.FontSize=16;
%theta4 position
Th4=uicontrol(f);
Th4.Style='edit';
Th4.Position=[80 160 50 25];
Th4.Units='pixels';
%theta4 label
LTH4=uicontrol(f);
LTH4.Style='text';
LTH4.Position=[30 160 40 25];
LTH4.String='TH4';
LTH4.FontSize=16;

    function Inverse(src,event)
        X=str2double(get(PX,'String'));
        Y=str2double(get(PY,'String'));
        Z=str2double(get(PZ,'String'));
        
        ik=IK(X,Y,Z);
        disp(ik)
        L(1)=Link([ 0  6.8 0    pi/2]);
        L(2)=Link([ 0  0   9.3  0   ]);
        L(3)=Link([ 0  0   10.9 0   ]);
        L(4)=Link([ 0  0   8.5  pi/2]);
        
        robot=SerialLink(L);
        
        
        figure(2);
        robot.name='A1';
        subplot(2,2,1);
        robot.plot(ik(1,1:4));
        robot.name='A2';
        subplot(2,2,2);
        robot.plot(ik(2,1:4));
        robot.name='A3';
        subplot(2,2,3);
        robot.plot(ik(3,1:4));
        robot.name='A4';
        subplot(2,2,4);
        robot.plot(ik(4,1:4));
        display(ik)
    end

    function Forward(src,event)
        t1=str2double(get(Th1,'String'));
        t2=str2double(get(Th2,'String'));
        t3=str2double(get(Th3,'String'));
        t4=str2double(get(Th4,'String'));
        
        L(1)=Link([ 0  6.8 0    pi/2]);
        L(2)=Link([ 0  0   9.3  0   ]);
        L(3)=Link([ 0  0   10.9 0   ]);
        L(4)=Link([ 0  0   8.5  pi/2]);
        
        robot=SerialLink(L);
        robot.name='AL5A';
        %dynamically ploting the robot
        if t1>=0
            for i=0:1:t1
                robot.plot([deg2rad(i) 0 0 0]);
            end
        end
        if t1<0
            T1=-1*t1;
            for i=0:1:T1
                robot.plot([deg2rad(-1*i) 0 0 0]);
            end
        end
        if t2>=0
            for j=0:1:t2
                robot.plot([deg2rad(t1) deg2rad(j) 0 0]);
            end
        end
        if t2<0
            T2=-1*t2;
            for j=0:1:T2
                robot.plot([deg2rad(t1) deg2rad(-1*j) 0 0]);
            end
        end
        if t3>=0
            for k=0:1:t3
                robot.plot([deg2rad(t1) deg2rad(t2) deg2rad(k) 0]);
            end
        end
        if t3<0
            T3=-1*t3;
            for k=0:1:T3
                robot.plot([deg2rad(t1) deg2rad(t2) deg2rad(-1*k) 0]);
            end
        end
        if t4>=0
            for s=0:1:t4
                robot.plot([deg2rad(t1) deg2rad(t2) deg2rad(t3) deg2rad(s)]);
            end
        end
        if t4<0
            T4=-1*t4;
            for s=0:1:T4
                robot.plot([deg2rad(t1) deg2rad(t2) deg2rad(t3) deg2rad(-1*s)]);
            end
        end
        
        table=uitable(f);
        table.Position=[320 5 300 160];
        table.Data=FK(deg2rad(t1) ,deg2rad(t2) ,deg2rad(t3) ,deg2rad(t4));
        
    end

end