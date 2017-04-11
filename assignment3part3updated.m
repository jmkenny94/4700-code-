%%Length of System
L = 200e-9

%%Width of System
W = 100e-9;

%%Number of Elements
Nx = 30
Ny = 20

%%Initial Voltage
Vo = 0.1

%%Discretize Xspace
xvec = linspace(0,L,Nx); %%start at 0, increment to L, length of vector is N
dx = xvec(2)-xvec(1);

%%Discretize Yspace
yvec = linspace(0,W,Ny);
dy = yvec(2)-yvec(1);

%%Allocate memory for voltage
V = zeros(length(xvec),length(yvec));

%% sigma plot
for idx = 1:Nx
    for idy = 1:Ny
        x = xvec(idx);
        y = yvec(idy);
        if x > 80e-9 & x < 120e-9 & y > -1 & y < 25e-9
            sigma(idx,idy) = 0.01;
        elseif x > 80e-9 & x < 120e-9 & y > 75e-9 & y < 110e-9
            sigma(idx,idy) = 0.01;
        else
            sigma(idx,idy) = 1;
        end
    end
end

R0 = 1./sigma;

%%Set up 4 resistor network for each point
for idx = 2:Nx-1
    for idy = 2:Ny-1
        R(idx,idy,1) = (R0(idx-1,idy)+R0(idx,idy))/2;
        R(idx,idy,2) = (R0(idx+1,idy)+R0(idx,idy))/2;
        R(idx,idy,3) = (R0(idx,idy+1)+R0(idx,idy))/2;
        R(idx,idy,4) = (R0(idx,idy-1)+R0(idx,idy))/2;
    end
end

%%initialize BCs
V(1,:) = Vo;


%%For loop to calculate V 
for n = 1:1000
    for idy = 1:Ny 
        for idx = 2:Nx-1 
            if idy == 1
                V(idx,idy) = V(idx,idy+1);
            elseif idy == Ny
                V(idx,idy) = V(idx,idy-1);
            else
                V(idx,idy) = ((((V(idx-1,idy))/R(idx,idy,1))+...
                    ((V(idx+1,idy))/R(idx,idy,2))+(V(idx,idy+1))...
                    /R(idx,idy,3))+((V(idx,idy-1))/R(idx,idy,4)))/...
                    ((1/R(idx,idy,1))+(1/R(idx,idy,2))+(1/R(idx,idy,3))+(1/R(idx,idy,4)));
            end
        end
    end
end

%E fields 
[ey,ex] = gradient(-V);

%define constants
m = 9.1e-31;
k = 1.381e-23;
T = 300;
np = 200000; %calculated to give density of 10^15 per cm squared 

xp = rand(np,1)*L;
yp = rand(np,1)*W;
%velocity vectors and std
vxstd = 0.5;
vx = sqrt(2*k*T/m)*randn(np,1)*vxstd;
vy = sqrt(2*k*T/m)*randn(np,1)*vxstd;

% %%Ex = 0.1V/(200*10^-9m)
% %%Ex = 500000V/m
% Ex = 500000

%%elementary charge
q = 1.602e-19

dt = 5e-15;



idx = zeros(np,1);
idy = zeros(np,1);
ax= zeros(np,1);
ay= zeros(np,1);

for t = 1:200
    tv(t) = t;
    
    idx = round(xp/dx)+1;
    idy = round(yp/dy)+1;
  
   
     for i=1:np
        b=idx(i);
        n=idy(i);
  
        
        %%Scaling factor here to convert efields from assignment to into
        %%V/nm. They were in V/m which was why acceleration wasnt enough to
        %%see curvature 
    ax(i) = 1.0e+9* q.*ex(b,n)/m;
     ay(i) = 1.0e+9* q.*ey(b,n)/m;
    end 
    
   
    dvx = ax*dt;
    
   
    dvy = ay*dt;
    
    VVx = (vx + dvx*t);
    VVy = (vy + dvy*t);
    
    x = xp + VVx*dt;
    y = yp + VVy*dt;
    
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
    
    for i = 1:np
        
        %LEFT WALL BOTTOM
        if (x(i)>75e-9) && (x(i)<85e-9) && (y(i)>0) && (y(i)<30e-9)
            vx(i) = -vx(i);
        end 
        
        %RIGHT WALL BOTTOM
        if (x(i)>115e-9) && (x(i)<125e-9) && (y(i)>0) && (y(i)<30e-9)
            vx(i) = -vx(i);
        end 
        
        %LOWER BOX CONNECT 
        if (x(i)>75e-9) && (x(i)<125e-9) && (y(i)>20e-9) && (y(i)<30e-9)
            vy(i) = -vy(i);
        end 
        
    end 
        
 
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
    pause(0.01)
    
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

 %logical indexing for density and temp plot 
    for ii = 1:20
        dx = L/20;
        %cell boundaries 
        x1 = (ii-1)*dx;
        x2 = x1 + dx;
       
        Ix = x > x1 & x < x2;
        
        for jj = 1:20
            dy = W/20;
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
    xlabel('Divisions of 10x10^-9 (meters)')
    ylabel('Divisions of 5x10^-9 (meters)')

