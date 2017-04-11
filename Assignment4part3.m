%%Vo with random noise 

R=20;
C=10*10^-6;
Vs=1;
Vsi = 1;
t=linspace(0,0.005,10000);
Vc=zeros(1,length(t));
RC=R*C;
tx = (t(2)-t(1));
Vc(length(t))=1;
I = (Vsi/R)*randn(1,length(t));
for i=1:length(t)-1
    Vc(i+1)=((C*(Vc(i))/tx)+(Vs/R))*(((C/tx)+1/R)^-1)*R;
    Vc(i+1) = Vc(i+1) + I(i+1);
    Vc(i+1) = Vc(i+1)/R;
    Vc2 = (Vc(i+1))^2
end
figure()
plot(t,Vc)
title('Output Voltage with Random Noise') 
xlabel('Time (s)')
ylabel('Voltage (V)')

V=fft(Vc)
figure()
plot(t,V)
title('fft Output ') 
xlabel('Time (s)')
ylabel('Voltage (V)')


VRMS = sum(Vc2(:))/((length(t)-1)-1);
VRMS = sqrt(VRMS);





