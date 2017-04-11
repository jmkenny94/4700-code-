clear variables
%define constants
m = 9.1e-31;
k = 1.381e-23;
T = 300;
np = 200000; %To give density of 10^15 per cm squared asked for in outline
%Width Length of system & starting positions
L = 200e-9;
W = 100e-9;
xp = rand(np,1)*L;
yp = rand(np,1)*W;
%velocity vectors and std
vxstd = 0.5;
vx = sqrt(2*k*T/m)*randn(np,1)*vxstd;
vy = sqrt(2*k*T/m)*randn(np,1)*vxstd;

%%Ex = 0.1V/(200*10^-9m)
%%Ex = 500000V/m
Ex = 500000;

%%elementary charge
q = 1.602e-19;

%%qE = ma --> a=qE/m
a = q*Ex/m;

dt = 5e-15;

dv = a*dt

for t = 1:200
    tv(t) = t;
    VV = (vx + dv*t);
    x = xp + VV*dt; 
    y = yp + vy*dt;
    %logical indexing of x boundaries
    ix = x > 200e-9;
    x(ix) = x(ix)-200e-9;
    ux = x < 0 ;
    x(ux) = x(ux) +200e-9;
    %logical indexing of x boundaries
    iy = y > 100e-9 ;
    vy(iy) = -vy(iy);
    uy = y < 0 ;
    vy(uy) = -vy(uy);
%     probability of scattering
    Pscat = (1-exp(-dt/0.2e-12));
    scat = rand(np,1);
    iscat = scat<Pscat;
    nones = nnz(iscat);
%     generate new velocity if scattered
    vx(iscat) = sqrt(2*k*T/m)*randn(nones,1)*vxstd;
    vy(iscat) = sqrt(2*k*T/m)*randn(nones,1)*vxstd;
    
    %Positive current contribution 
    Ipos = nnz(ix);
    
    %%Negative contribution to current 
    Ineg = nnz(ux);
    
    %%Current Calculation 
    current(t) = (Ipos-Ineg)*q/dt;
    time(t)= t;
    
    figure(1)
    plot(x(1:20000:end),y(1:20000:end),'o');
    grid on
    hold on
    title('Electron Mapping')
    xlabel('x Position, x (m)')
    ylabel('Y Position, y (m)')
    axis ([0 L 0 W])
    
    
    msv = mean(vx.^2 + vy.^2);
    temp(t) = (m*msv/(k));
    %figure(2)
    %plot(tv,temp)
    %title('Temperature vs. Time')
    %xlabel('Time')
    %ylabel('Temperature(C)')
    
    
    xp = x;
    yp = y;
    
    pause(0.1)
end

figure (2)
plot(time,current)
title('Current Vs. Time')
ylabel('Current (A)')
xlabel('Time (s)')

    %logical indexing for density and temp plot 
    for ii = 1:10
        dx = L/10;
        %cell boundaries 
        x1 = (ii-1)*dx;
        x2 = x1 + dx;
       
        Ix = x > x1 & x < x2;
        
        for jj = 1:5
            dy = W/5;
            %cell boundaries
            y1 = (jj-1)*dy;
            y2 = y1 + dy;
            
            Iy = y > y1 & y < y2;
            %particles in cell 
            Ixy = Ix&Iy;
            %positions
            xpa = x(Ixy);
            ypa = y(Ixy);
            %velocities
            vxpa = vx(Ixy);
            vypa = vy(Ixy);
            % number of particles in cell
            npa = length(xpa);
            
            if npa == 0
                Temp(ii,jj) = 300;
                %density plot 
                elecdense(ii,jj)= npa./(W*L);
            else
                %temp plot 
                Temp(ii,jj) = m*sum(vxpa.^2 + vypa.^2)/(2*k*npa);
                elecdense(ii,jj)= npa./(W*L);
                
            end
        end 
    end 
    
    figure(4)
    pcolor(Temp)
    title('Surface Temperature Plot')
    zlabel('Temperature (C)')
    xlabel('Divisions of 20x10^-9 (meters)')
    ylabel('Divisions of 20x10^-9 (meters)')  
    
    figure (5)
    pcolor(elecdense)
    title('Electron Density Map')
    zlabel(' e- Density (e-/m^2)')
    xlabel('Divisions of 20x10^-9 (meters)')
    ylabel('Divisions of 20x10^-9 (meters)')
    
    
