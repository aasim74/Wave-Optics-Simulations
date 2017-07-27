%Code to simulate X-ray diffraction by first animating the incoming wave
%before reaching the aperture and then animating the wave after it passes
%through the aperture and hitting the screen

%Define the variables
c = 3e8;
lambda = 630e-9;
f = c/lambda;
k = 2*pi/lambda;
omega = 2*pi*f;
A = 1;
t=0;

x = linspace(-15e-6, 15e-6, 501);
y = linspace(-15e-6, 15e-6, 501);

incoming_wave = zeros(length(x), length(y));

size(incoming_wave)

%Creates a 2-D matrix of intensities at different x and y values
for i = 1: 1: length(x)
    for j = 1: 1: length(y)
        incoming_wave(i,j) = A*exp(1j*k*(x(i)+y(j)) - omega*t);
    end
end

size(incoming_wave)

t = 0: 1: 100;

%Creates a 3-D array with time being the third dimension
for l = 1: 1: length(t)
    incoming_wave(:, :, l) = A*exp(1j*k*(x(i)+y(j)) - omega*t(l));
end

figure(1);

%Plot the wave as a function of x-y position and evolving with time t
for i = 1: 1: length(x)
    for j = 1: 1: length(y)
        for t=0
            surf(x, y, real(incoming_wave))
        end
    end
end


size(incoming_wave)