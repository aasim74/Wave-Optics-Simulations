%FDTD Code Brian Egenriether
%This code is a finite-difference time-domain simulation of 3 gaussian
%pulses in a resonant cavity with perfectly conducting walls (PEC)
%It creates an animation of the results but if you want to let the system
%get to steady state you'll need around 50,000 steps (i.e.nMax= 50000)
%This cannot be animated, you will run out of memory. Comment out the
%getframe command at the end of the outer loop and uncomment the plots at
%the bottom to get slices of the steady state result and the FFT.

clc
clear
nMax = 700;  % number of time steps
scalefactor=1;  %Granularity (for animations don't exceed 2)


mu =  1.2566306E-6;
eps0 = 8.854E-12;
epsr = 4.2;    
eps = eps0*epsr;
c = 1/sqrt(mu*eps);
cfln = 0.995; %stability criterion
tw = 15;% %half width of pulse  10-20 looks nice
n0 = 5*tw; %time delay for sources (stability criterion)
Lx = 0.08;  %dimensions of resonant cavity
Ly = 0.06;
Lz = 0.048;

nx = 80*scalefactor;  %3 equal discretizations
ny =60*scalefactor;
nz= 48*scalefactor;

dx = Lx/nx;   %finite differences for differentials
dy = Ly/ny;
dz = Lz/nz;
dt = cfln/(c*sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2));%time stability criterion 
fres = c/sqrt(Lx^2+Ly^2+Lz^2);
%FDTD Coefficients
cHxy = dt/(mu*dy);
cHxz = dt/(mu*dz);
cHyz = dt/(mu*dz);
cHyx = dt/(mu*dx);
cHzx = dt/(mu*dx);
cHzy = dt/(mu*dy);
cExy = dt/(eps*dy);
cExz = dt/(eps*dz);
cEyz = dt/(eps*dz);
cEyx = dt/(eps*dx);
cEzx = dt/(eps*dx);
cEzy = dt/(eps*dy);

%PEC Boundary conditoins
%Set all tangential E fields to zero and never update them.  This will keep
%them zero. Update loops only modify inner cells
Hx=zeros(nx,ny,nz);
Hy=zeros(nx,ny,nz);
Hz=zeros(nx,ny,nz);
Ex=zeros(nx,ny,nz);
Ey=zeros(nx,ny,nz);
Ez=zeros(nx,ny,nz);

%=====================Outer Loop
for n=1:nMax
   
 

%==================== Main update loops for the Electric fields:
%===================Ex
  for k = 2:nz-1
      for j = 2:ny-1
          for i = 1:nx-1 %1->2
              Ex(i,j,k) = Ex(i,j,k)+cExy*(Hz(i,j,k)-Hz(i,j-1,k))-cExz*(Hy(i,j,k)-Hy(i,j,k-1));                                  
          end
      end
  end
  %===================Ey
  for k = 2:nz-1
      for j = 1:ny-1%1->2
          for i = 2:nx-1
              Ey(i,j,k) = Ey(i,j,k)+cEyz*(Hx(i,j,k)-Hx(i,j,k-1))-cEyx*(Hz(i,j,k)-Hz(i-1,j,k));
          end
      end
  end  
  %===================Ez
  for k = 1:nz-1  %1->2
    for j = 2:ny-1
        for i = 2:nx-1
            Ez(i,j,k) = Ez(i,j,k)+cEzx*(Hy(i,j,k)-Hy(i-1,j,k))-cEzy*(Hx(i,j,k)-Hx(i,j-1,k));
        end
    end
  end

%Define E field components at 3 seperate arbitrary points with Gaussian pulses 
Ex(30*scalefactor,20*scalefactor,25*scalefactor)=2*tw*exp(-((n-n0)/tw)^2); 
Ey(25*scalefactor,45*scalefactor,20*scalefactor)=2*tw*exp(-((n-n0)/tw)^2);  
Ez(55*scalefactor,30*scalefactor,23*scalefactor)=2*tw*exp(-((n-n0)/tw)^2); 

%==================== Main update loops for the Magnetic fields:
%===================Hx
for k = 1:nz-1  %1->2 all 3
    for j = 1:ny-1
        for i = 1:nx %-1
            Hx(i,j,k) = Hx(i,j,k)-cHxy*(Ez(i,j+1,k)-Ez(i,j,k))+cHxz*(Ey(i,j,k+1)-Ey(i,j,k));                                  
        end
    end
end
%===================Hy
for k = 1:nz-1  %1->2 all 3
    for j = 1:ny %-1
        for i = 1:nx-1
            Hy(i,j,k) = Hy(i,j,k)-cHyz*(Ex(i,j,k+1)-Ex(i,j,k))+cHyx*(Ez(i+1,j,k)-Ez(i,j,k));            
        end
    end
end
%===================Hz
for k = 1:nz %-1  %1->2 all 3
    for j = 1:ny-1
        for i = 1:nx-1
             Hz(i,j,k) = Hz(i,j,k)-cHzx*(Ey(i+1,j,k)-Ey(i,j,k))+cHzy*(Ex(i,j+1,k)-Ex(i,j,k));   
        end
    end
end    


%Define H field components at 3 seperate arbitrary points with Gaussian pulse  
%Must put this after the H field update for stability
% Hx(30,30,30)=exp(-((n-n0)/tw)^2);
% Hy(45,35,20)=exp(-((n-n0)/tw)^2);  
% Hz(25,40,20)=exp(-((n-n0)/tw)^2);  
clc
round(n/nMax*100)   %gives progress update percentage on home screen

surf(Ez(:,:,28*scalefactor)); 
set(gcf,'renderer','zbuffer')
%hold on
whitebg('black');
grid off
set(gcf,'Position',[20 50 1500 800]);
axis([-10*scalefactor 60*scalefactor -20*scalefactor 70*scalefactor -.01/scalefactor .01/scalefactor])
 view([3,4,5])
 
Mov(n)=getframe; % comment this out for nonanimated analysis
   
   
end
movie(Mov,20,40);    % comment this out for nonanimated analysis