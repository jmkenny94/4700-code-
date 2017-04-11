%Part_B
R=20;
C=10*10^-6;
E=1;
tt=linspace(0,0.001,100);
Vc=zeros(100,1);

for i=1:100
syms t
I=int(sin(10^6*t)/R*(exp(-t/(R*C))),t,[0,((6.283*10^-6*30)/100)*i]);
Vc(i)=1/C*(I);
end 


plot(tt,Vc)
V=fft(Vc)
title('Output when Using Sin Wave Input') 
xlabel('Time (s)')
ylabel('Voltage (V)')



figure()
plot(tt,V)
title('fft Output ') 
xlabel('Time (s)')
ylabel('Voltage (V)')
