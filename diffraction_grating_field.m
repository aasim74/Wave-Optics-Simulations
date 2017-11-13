%{
Fresnel Diffraction through a rectangular aperture. Solutions were taken
out from the Goodmans Optics Textbook and transmitted Intensity plotted.

Aasim Patel
%}

%Define the constants
n=10;
D = 10e-6; 
wavelength = 630e-9;
delta_z = 10;
k = 2*pi/wavelength;

xi = 0:0.1:10;
yi = 0:0.1:10;
U = zeros(length(xi), length(yi));

tic;
startTime = tic;
        
        alpha1 = -sqrt(2/wavelength*delta_z)*(D/2 + xi);
        alpha2 = sqrt(2/wavelength*delta_z)*(D/2 - xi);
        beta1 = -sqrt(2/wavelength*delta_z)*(D/2 + yi);
        beta2 = sqrt(2/wavelength*delta_z)*(D/2 - yi);
        
        C_alpha1 = fresnelc(alpha1);
        C_alpha2 = fresnelc(alpha2);
        C_beta1  = fresnelc(beta1);
        C_beta2 = fresnelc(beta2);
        
        S_alpha1 = fresnels(alpha1);
        S_alpha2 = fresnels(alpha2);
        S_beta1  = fresnels(beta1);
        S_beta2 = fresnels(beta2);
        
for i = 1: length(xi)
    for j = 1: length(yi)
        
        U(i,j) = (exp(1i*k*delta_z)/2*1i) * ((C_alpha2(i) - C_alpha1(i)).^2 + 1i*(S_alpha2(i) - S_alpha1(i)).^2)...
        .*((C_beta2(i) - C_beta1(i)).^2 + 1i*(S_beta2(i) - S_beta1(i)).^2);
    end
end

figure(1);
plot(xi, fft(real(U)))
title('Transmitted Field');
xlabel('x');
ylabel('\psi_{t}');

figure(2);
plot(xi, fft(U*conj(U)));
title('Transmitted Intensity');
xlabel('x');
ylabel('I(x)');

endTime = toc();

timeElapsed = endTime - startTime;
disp(timeElapsed);