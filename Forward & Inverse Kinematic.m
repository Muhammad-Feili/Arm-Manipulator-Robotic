function RoboticPorject
clear all
close all
clc


    function dh=A(th, d, a, alpha)
             dh=[cos(th),   -sin(th)*cos(alpha),  sin(th)*sin(alpha),   a*cos(th);  
                sin(th),   cos(th)*cos(alpha),   -cos(th)*sin(alpha),  a*sin(th);
                0,         sin(alpha),           cos(alpha),           d        ;
                0,         0,                    0,                    1       ];
    end




    function Hdir=ForwardKinematics(th1,th2,th3,th4,th5,th6)
        
        C1=cos(th1); C2=cos(th2); C4=cos(th4); C5=cos(th5); C6=cos(th6);
        S1=sin(th1); S2=sin(th2); S4=sin(th4); S5=sin(th5); S6=sin(th6);
        C23=cos(th2+th3); S23=sin(th2+th3);
        
        a1=0; a2=0.7; d3=0.0948;
        
        nx=C1*(C23*(C4*C5*C6-S4*S6)-S23*S5*C6)+S1*(S4*C5*C6+C4*S6);
        ny=S1*(C23*(C4*C5*C6-S4*S6)-S23*S5*C6)-C1*(S4*C5*C6+C4*S6);
        nz=S23*(C4*C5*C6-S4*S6)+C23*S5*C6;
        ox=C1*(-C23*(C4*C5*C6+S4*C6)+S23*S5*S6)+S1*(C4*C6-S4*S5*S6);
        oy=S1*(-C23*(C4*C5*C6+S4*C6)+S23*S5*S6)-C1*(C4*C6-S4*S5*S6);
        oz=-S23*(C4*C5*C6+S4*C6)-C23*S5*S6;
        ax=C1*(C23*C4*S5+S23*C5)-C1*S4*S5;
        ay=S1*(C23*C4*S5+S23*C5)-C1*S4*S5;
        az=S23*C4*S5-C23*C5;
        px=C1*(C2*a2+a1)+S1*d3;
        py=S1*(C2*a2+a1)-C1*d3;
        pz=S2*a2;
        
        Hdir=[ nx ox ax px; ny oy ay py; nz oz az pz; 0 0 0 1];
        
    end
    function Hinv=InverseKinematics(Hdir,Pos)
        
        nx=Hdir(1,1); ox=Hdir(1,2); ax=Hdir(1,3); px=Pos(1,1); %Hdir(1,4);
        ny=Hdir(2,1); oy=Hdir(2,2); ay=Hdir(2,3); py=Pos(1,2); %Hdir(2,4);
        nz=Hdir(3,1); oz=Hdir(3,2); az=Hdir(3,3); pz=Pos(1,3); %Hdir(3,4);
        
        a2=0.7; d3=0.0948; a1=0;
        
        th2=asin(pz/a2); C2=cos(th2);
        th1=acos((px*(a2*C2+a1)-py*d3)/(px^2+py^2)); C1=cos(th1); S1=sin(th1);
        th23=acos(sqrt((pz^2-a2^2+(px*C1+py*S1-a1)^2)/(2*(pz^2)))); C23=cos(th23); S23=sin(th23);
        th3=th23-th2;
        th4=atan((ax*S1-ay*C1)/(ax*C1*C23+ay*S1*S23+az*S23)); C4=cos(th4); S4=sin(th4);
        th5=-atan((az*C1*C23*C4+ay*S1*S23*C4+az*S1*S4-ay*C1*S4)/(az*C23-ay*S1*S23-az*S23));
        th6=atan((oz*C23-ox*C1*S23-oy*S1*S23)/(nx*C1*S23+ny*S1*S23-nz*C23));
        
        Hinv=[th1 th2 th3 th4 th5 th6];
        
        
    end
% values for simulataing
%6 DOF robot manipulator
entVals=[526.628 139.692 1100.599;
        527.0173 140.1437 1100.13;
        529.6739 143.2708 1096.88;
        536.4599 151.6324 1088.162;
        548.2687 167.6693 1071.309;
        564.3856 193.6511 1043.588;
        581.9591 231.3008 1002.315;
        595.9011 280.9468 945.3523;
        599.6258 340.4005 871.9466;
        586.8507 404.2318 783.6389;
        554.0557 464.3222 684.7641;
        502.4461 512.0784 582.0885;
        438.1163 541.5359 483.4645;
        370.0438 551.6382 395.9441;
        306.9984 546.2927 324.1966;
        255.0976 532.2926 269.9224;
        216.9997 516.4974 232.3015;
        192.4576 503.7352 208.942;
        179.3259 496.0306 196.7058;
        174.3527 492.9474 192.118]*10^(-3);
    
    
    
jointVals=[0        0         0       0        0        0;
           0.049105 0.049105 -0.04911 0.049105 0.049105 0.049105;
           0.387064 0.387064 -0.38706 0.387064 0.387064 0.387064;
           1.274457 1.274457 -1.27446 1.274457 1.274457 1.274457;
           2.918079 2.918079 -2.91808 2.918079 2.918079 2.918079;
           5.450703 5.450703 -5.4507  5.450703 5.450703 5.450703;
           8.918079 8.918079 -8.91808 8.918079 8.918079 8.918079;
           13.27446 13.27446 -13.2745 13.27446 13.27446 13.27446;
           18.38706 18.38706 -18.3871 18.38706 18.38706 18.38706;
           24.04911 24.04911 -24.0491 24.04911 24.04911 24.04911;
           30       30       -30      30       30       30;
           35.95089 35.95089 -35.9509 35.95089 35.95089 35.95089;
           41.61294 41.61294 -41.6129 41.61294 41.61294 41.61294;
           46.72554 46.72554 -46.7255 46.72554 46.72554 46.72554;
           51.08192 51.08192 -51.0819 51.08192 51.08192 51.08192;
           54.5493  54.5493  -54.5493 54.5493  54.5493  54.5493;
           57.08192 57.08192 -57.0819 57.08192 57.08192 57.08192;
           58.72554 58.72554 -58.7255 58.72554 58.72554 58.72554;
           59.61294 59.61294 -59.6129 59.61294 59.61294 59.61294;
           59.95089 59.95089 -59.9509 59.95089 59.95089 59.95089]; 

%intial values of variable joints
%ForwardKinematics(th1,th2,th3,th4,th5,th6)


%TH1=deg2rad(-185); TH2=deg2rad(-155); TH3=deg2rad(-130);
%TH4=deg2rad(-350); TH5=deg2rad(-130); TH6=deg2rad(-350);

%hDIR=ForwardKinematics(TH1,TH2,TH3,TH4,TH5,TH6);
%hDIR=ForwardKinematics(-185,-155,-130,-350,-130,-350);
 
simVals=zeros(20,3);

for i =1:20
    B=ForwardKinematics(deg2rad(jointVals(i,1)),deg2rad(jointVals(i,2)),deg2rad(jointVals(i,3)),deg2rad(jointVals(i,4)),deg2rad(jointVals(i,5)),deg2rad(jointVals(i,6)));
    simVals(i,1:3)=B(1:3,4);
end

X=simVals(1:20,1);
Y=simVals(1:20,2);
Z=simVals(1:20,3);

disp(simVals)

plot(jointVals(1:20,1),X)
hold on
plot(jointVals(1:20,1),Y)
hold on
plot(jointVals(1:20,1),Z)


end