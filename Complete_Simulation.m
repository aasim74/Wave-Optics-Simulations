%Code to simulate X-ray diffraction by first animating the incoming wave
%before reaching the aperture and then animating the wave after it passes
%through the aperture and hitting the screen

%Define the variables------------------------------------------------------
c = 3e8;
lambda = 630e-9;
f = c/lambda;
k = 2*pi/lambda;
omega = 2*pi*f;
A = 1;

x = linspace(-1e-6, 1e-6, 201);
y = linspace(-1e-6, 1e-6, 201);
t= 0: 1: 10;


%Creates a 2-D matrix of intensities at different x and y values-----------
incoming_wave = zeros(length(x), length(y));

for i = 1: 1: length(x)
    for j = 1: 1: length(y)
        for l = 1: 1: length(t)
        incoming_wave(i, j, l) = A*exp(1j*k*(x(i)+y(j) - omega*t(l)));
        end
    end
end

%Plot 2-D graph------------------------------------------------------------
figure(1);
surf(x, y, real(incoming_wave(:,:,1)))
xlabel('x-axis');
ylabel('y-axis');
zlabel('Amplitude');