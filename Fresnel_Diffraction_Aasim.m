%{
--------------------------------------------------------------------------
Program to find the diffraction pattern using the Fresnel Approximation
The program uses Matlabs built in functions to solve the fresnel integrals
to solve the C(x) integral and S(x) integrals which were used to find the 
intensity
--------------------------------------------------------------------------
Written by Aasim Patel
26 July 2017
--------------------------------------------------------------------------
%}

%Define the aperture -----------------------------------------------------
W = 10e-6;
xi = 0: 10e-7 :10e-7;   % xi-axis
eta = 0: 10e-7: 10e-5;  % eta-axis
lambda = 630e-9;
z = 1; % distance from the aperture
NF = W^2/lambda*z;
%prompt = 'Please input the Fresnel Number NF: ';
%NF = input(prompt) ;% Fresnel Number
x = linspace(-15e-6,15e-6,501); %observation region
y = linspace(-15e-6,15e-6,501);
X = x/lambda*z; % Normalised X variable
Y = y/lambda*z; % Normalised Y variable

StartClock = tic;

%find the alpha and beta values ------------------------------------------
alpha1 = -sqrt(2)*(sqrt(NF) + X);
alpha2 = -sqrt(2)*(sqrt(NF) - X);
beta1 = -sqrt(2)*(sqrt(NF) + Y);
beta2 = -sqrt(2)*(sqrt(NF) - Y);

%find the fresnel integral values ----------------------------------------
C_alpha1 = fresnelc(alpha1);
C_alpha2 = fresnelc(alpha2);
C_beta1 = fresnelc(beta1);
C_beta2 = fresnelc(beta2);
S_alpha1 = fresnels(alpha1);
S_alpha2 = fresnels(alpha2);
S_beta1 = fresnels(beta1);
S_beta2 = fresnels(beta2);

%find the relative intensities -------------------------------------------
intensity = 0.25*(((C_alpha2 - C_alpha1).^2 +1j*(S_alpha2 - S_alpha1).^2).*...
    ((C_beta2 - C_beta1).^2 + 1j*(S_beta2 - S_beta1).^2));

%plot --------------------------------------------------------------------
figure(1);
plot(x, real(intensity))
xlabel('x-displacement on screen')
ylabel('Relative Intensity at fixed distance z')
title('Intensity as a function of prosition in x')

figure(2);
plot(intensity, y)
ylabel('y-displacement on screen')
xlabel('Relative Intensity at fixed distance z')
title('Intensity as a function of y')
%title('Intensity as a function of position in y')

intensity = zeros(length(x), length(y));

for i = 1: 1: length(x)
    for j = 1: 1: length(y)
        
    intensity(i,j) =  0.25 * ( ( (C_alpha2(i) - C_alpha1(i) ).^2 ...
    + 1j*( S_alpha2(i) - S_alpha1(i) ).^2) .* ...
    (( C_beta2(j) - C_beta1(j) ).^2 + 1j*( S_beta2(j) - S_beta1(j) ).^2));

    end
end

figure(3);
imshow(real(intensity))
StopClock = toc(StartClock);
disp(StopClock)

% End --------------------------------------------------------------------