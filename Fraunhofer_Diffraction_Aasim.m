%{
Program to find the Fraunhofer Diffraction Pattern for a range of x and
y values for a given width and length of a rectangular/square aperture. 

Written by Aasim Patel
27 July 2017

%}

% Define values------------------------------------------------------------
lambda = 630e-9;    % Wavelength
z = 100;    % Distance
Wx = 1e-2;  % Aperture Width x-axis
Wy = 1e-2;  % Aperture Length y-axis
A = 4*Wx*Wy; % Aperture Area

x = linspace(-15e-7, 15e-7, 501);
y = linspace(-15e-7, 15e-7, 501);

intensity = (A^2/lambda*z) * sinc(2*Wx*x/lambda*z).^2 ...
            .* sinc(2*Wy*y/lambda*z).^2;

%Plot----------------------------------------------------------------------
figure(1);
plot(x, intensity)
title('Intensity in the x-axis');
xlabel('x-displacement');
ylabel('Relative Intensity');

figure(2);
plot(intensity, y)
title('Intensity in the y-axis');
xlabel('Relative Intensity');
ylabel('y-displacement');

intensity = zeros(length(x), length(y));

for i = 1 : 1: length(x)
    for j = 1: 1: length(y)
        intensity(i,j) = (A^2/lambda*z) * sinc(2*Wx*x(i)/lambda*z).^2 ...
        .* sinc(2*Wy*y(j)/lambda*z).^2;
    end
end

figure(3);
imshow(intensity)
title('Fraunhofer Diffraction Pattern');

%End-----------------------------------------------------------------------