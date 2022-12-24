clc
clear

%=================%
%Identity%
%=================%
% Name   : Muhammad Fadli Alfarizi
% NIM    : 13121140

%=================%
%Input%
%=================%
XYZ = "140"; % NIM

%=================%
%System Parameters%
%=================%
strXYZ = char(XYZ);
X = strXYZ(1);
Y = strXYZ(2);
Z = strXYZ(3);

mod1 = mod(str2num(XYZ),3);
if mod1 == 0
    H = 22.5*(10^-2); % in meter
elseif mod1 == 1
    H = 25*(10^-2); % in meter
else
    H = 25*(10^-2); % in meter
end

mod2 = mod(str2num(XYZ),2);
% angular velocity
if mod2 == 0
    w2 = 4; % in rad/s
else
    w2 = 2; % in rad/s
end

L2=(30+str2num(Z))*(10^-2) ; % in meter
L3=(136+str2num(X))*(10^-2); % in meter
L5=(57 + str2num(Y)/4)*(10^-2); % in meter

% point O
uPjO1=[0;0];
uPjO2=[-L2/2;0];

% point B
uPjB2=[L2/2;0];
uPjB3=[-L3/2;0];

% point C
uPjC3=[-L3/6; 0];
uPjC5=[-L5/2; 0];

% point D
uPjD3=[L3/2;0];
uPjD4=[0;0];

% point F
uPjF5=[L5/2;0];
uPjF6=[0;0];

%==============%   
%Initial Values%
%==============%`
q=zeros(18,1);
q6init=pi/6;
%==============%

t=0:0.001:2;

NumData=size(t,2);
q_alltime=zeros(18,NumData);
v_alltime=zeros(18,NumData);
a_alltime=zeros(18,NumData);

for j=1:NumData
 
    epsilon=1e-12;
    delta_q_norm=1;
    i_newton=1;

    while abs(delta_q_norm)>epsilon

    %Transformation Matrix & Its derivative%
    A1=[cos(q(3)) -sin(q(3)); sin(q(3)) cos(q(3))]; % link 1
    dA1t1=[-sin(q(3)) -cos(q(3)); cos(q(3)) -sin(q(3))]; % link 1
    A2=[cos(q(6)) -sin(q(6)); sin(q(6)) cos(q(6))]; % link 2
    dA2t2=[-sin(q(6)) -cos(q(6)); cos(q(6)) -sin(q(6))]; % link 2
    A3=[cos(q(9)) -sin(q(9)); sin(q(9)) cos(q(9))]; % link 3
    dA3t3=[-sin(q(9)) -cos(q(9)); cos(q(9)) -sin(q(9))]; % link 3
    A4=[cos(q(12)) -sin(q(12)); sin(q(12)) cos(q(12))]; % link 4
    dA4t4=[-sin(q(12)) -cos(q(12)); cos(q(12)) -sin(q(12))]; % link 4
    A5=[cos(q(15)) -sin(q(15)); sin(q(15)) cos(q(15))]; % link 5
    dA5t5=[-sin(q(15)) -cos(q(15)); cos(q(15)) -sin(q(15))]; % link 5
    A6=[cos(q(18)) -sin(q(18)); sin(q(18)) cos(q(18))]; % link 6
    dA6t6=[-sin(q(18)) -cos(q(18)); cos(q(18)) -sin(q(18))]; % link 6

    %Constraint Equation%
    C=zeros(18,1);
    C45=[q(1);q(2)]+A1*uPjO1-[q(4);q(5)]-A2*uPjO2; % point O, Constrain Equation Matrix row 4-5
    C67=[q(4);q(5)]+A2*uPjB2-[q(7);q(8)]-A3*uPjB3; % point B, Constrain Equation Matrix row 6-7
    C89=[q(7);q(8)]+A3*uPjD3-[q(10);q(11)]-A4*uPjD4; % point D, Constrain Equation Matrix row 8-9
    C1011=[q(7);q(8)]+A3*uPjC3-[q(13);q(14)]-A5*uPjC5; % point C, Constrain Equation Matrix row 10-11
    C1213=[q(13);q(14)]+A5*uPjF5-[q(16);q(17)]-A6*uPjF6; % point F, Constrain Equation Matrix row 12-13
    C1415=[q(11);q(12)]; % slider 4, Constrain Equation Matrix row 14-15
    C1617=[q(17)-H;q(18)]; %slider 6, Constrain Equation Matrix row 16-17
    C12=q(6)-q6init-w2*t(j); % link 2, Constrain Equation Matrix row 18
    
    %Constrain Equation Matrix%
    C=[q(1);q(2);q(3);C45;C67;C89;C1011;C1213;C1415;C1617;C12];
    
    C_norm=sqrt(sum(C(1:18).^2))/18;

    %Constrain Jacobian Matrix%
    unit22=[1 0; 0 1];
    
    %Zeros 2x2 Matrix
    zeros22=[0 0; 0 0];

    %Ground Constraint%
    C0q=[1 zeros(1,17); 0 1 zeros(1,16); 0 0 1 zeros(1,15)]; % Constrain Jacobian Matrix row 1-3

    %Pin Join%
    CPjO=[zeros22 dA1t1*uPjO1 -unit22 -dA2t2*uPjO2 zeros(2,12)]; % Constrain Jacobian Matrix row 4-5
    CPjB=[zeros(2,3) unit22 dA2t2*uPjB2 -unit22 -dA3t3*uPjB3 zeros(2,9)]; % Constrain Jacobian Matrix row 6-7
    CPjD=[zeros(2,6) unit22 dA3t3*uPjD3 -unit22 -dA4t4*uPjD4 zeros(2,6)]; % Constrain Jacobian Matrix row 8-9
    CPjC=[zeros(2,6) unit22 dA3t3*uPjC3 zeros(2,3) -unit22 -dA5t5*uPjC5 zeros(2,3)]; % Constrain Jacobian Matrix row 10-11
    CPjF=[zeros(2,12) unit22 dA5t5*uPjF5 -unit22 -dA6t6*uPjF6]; % Constrain Jacobian Matrix row 12-13

    %Slider%
    Cslider4=[zeros(1,10) 1 zeros(1,7); zeros(1,11) 1 zeros(1,6)]; % Constrain Jacobian Matrix row 14-15
    Cslider6=[zeros(1,16) 1 zeros(1,1); zeros(1,17) 1]; % % Constrain Jacobian Matrix row 16-17

    %Driving Constraint, Link 2% 
    Cdrive=[zeros(1,5) 1 zeros(1,12)]; % Constrain Jacobian Matrix row 18
 
    % Final Constrain Jacobian Matrix
    Cq=[C0q;CPjO;CPjB;CPjD;CPjC;CPjF;Cslider4;Cslider6;Cdrive];
   
    delta_q=inv(Cq)*(-C);
    delta_q_norm=sqrt(sum(delta_q(1:12).^2))/12;
    
    q=q+delta_q;

    i_newton=i_newton+1;

      if  i_newton > 30 %limit newton rhapson iteration in case non-convergence%
        break
      end  
    end

    Ct=[zeros(17,1);-w2];
    velocity=inv(Cq)*(-Ct);
    
    QdPjO=A1*uPjO1*(velocity(3,1))^2-A2*uPjO2*(velocity(6,1))^2;
    QdPjB=A2*uPjB2*(velocity(6,1))^2-A3*uPjB3*(velocity(9,1))^2;
    QdPjD=A3*uPjD3*(velocity(9,1))^2-A4*uPjD4*(velocity(12,1))^2;
    QdPjC=A3*uPjC3*(velocity(9,1))^2-A5*uPjC5*(velocity(15,1))^2;
    QdPjF=A5*uPjF5*(velocity(15,1))^2-A6*uPjF6*(velocity(18,1))^2;
    Qd=[zeros(3,1);QdPjO;QdPjB;QdPjD;QdPjC;QdPjF;zeros(5,1)];

    accel=inv(Cq)*Qd;

    q_alltime(:,j)=q;
    v_alltime(:,j)=velocity;
    a_alltime(:,j)=accel;

end

C_norm;
delta_q_norm;
i_newton;


figure(1);
plot(t,q_alltime(16,:));
grid on;
xlabel('Time [s]');
ylabel('Slider Position [m]');
title('Position of Slider 6 vs Time');

figure(2);
plot(t,v_alltime(16,:));
grid on;
xlabel('Time [s]');
ylabel('Slider Velocity [m/s]');
title('Velocity of Slider 6 vs Time');

figure(3);
plot(t,a_alltime(16,:));
grid on;
xlabel('Time [s]');
ylabel('Slider Velocity [m/s]');
title('Acceleration of Slider 6 vs Time');

figure(4);
plot(t,v_alltime(15,:));
grid on;
xlabel('Time [s]');
ylabel('Angular Velocity [m/s]');
title('Angular Velocity of Link 5 vs Time');

figure(5);
plot(t,a_alltime(15,:));
grid on;
xlabel('Time [s]');
ylabel('Angular Acceleration [m/s]');
title('Angular Acceleration of Link 5 vs Time');
