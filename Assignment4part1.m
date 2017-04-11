% Define Voltage source, Capacitor, and Time lapse
Vs = 1;
C = 10e-6
t = 0: 0.00001: 0.0010;

% Define resistor, time constant, and output
R = 20;
tau = R*C;
Vo = Vs * (1 - exp(-t/tau));

%Plot Transient Analysis
figure(1)
plot(t,Vo)
grid on
title('Transient Analysis RC circuit')
xlabel('Time(s)')
ylabel('Voltage Across Cap (V)')




%Transfer Function

f = 0:10:10000
w = 2*pi*f
H = 1./((1i*w*R*C)+1);

figure(2)
plot (f,abs(H))

figure(3)
plot(log(f),20*log10(H))
title('Transfer Function Frequency Response')
xlabel('Log Frequency (Hz)')
ylabel('dB Output')


%%Numerical Solving 
R=20;
C=10*10^-6;
Vs=1;
t=linspace(0,0.001,10000);
Vc=zeros(1,length(t));
RC=R*C;
tx = (t(2)-t(1));
Vc(length(t))=1;
for i=2:length(t)-2
    Vc(i+1)=((C*(Vc(i))/tx)+(Vs/R))*(((C/tx)+1/R)^-1);
    
end
figure()
plot(t,Vc)
title('Output From Numerical Method')
xlabel('Time(s)')
ylabel('Voltage Across Cap (V)')





%%Vrms Calculation --> Vrms = sqrt(4kTRB)

%Vrms = 0.000000045 V


