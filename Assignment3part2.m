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


%%Sigma Plot 
figure(1)
surf(sigma)
title('Sigma Plot')
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Sigma Value (1/R)')

%%Resistance Plot 
figure(2)
surf(R0)
title('Resistance Plot')
xlabel('X position (m)')
ylabel('Y position (m)')
zlabel('Resistance (ohms)')


%E fields 
[ey,ex] = gradient(-V);

[xx,yy] = meshgrid(yvec,xvec);

%%Electric Fields Plot 
figure (3)
contour(xx,yy,V)
hold on
quiver(yvec,xvec,ey,ex)
hold off
title('Electric Field Plot (Gradient(-V))(V/m)')
xlabel('Y position (m)')
ylabel ('X position (m)')


figure (4)
surface(V)
title('Potential Surface Plot')
xlabel('X Direction')
ylabel('Y Direction')
zlabel('Potential (V)')
