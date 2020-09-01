%% INITIALIZE MATLAB
clear all
clc
close all
format long

%% GRINGARTEN TYPE CURVE

% INCREMENTS
% time saving increments
tD2 = 0.1:0.01:0.99;
tD3 = 1:.1:9.9;
tD4 = 10:1:99;
tD5 = 100:10:999;
tD6 = 1000:100:1e4;

% tD2 = 0.1:0.01:0.99;
% tD3 = 1:.1:9.9;
% tD4 = 10:.1:99.9;
% tD5 = 100:1:999;
% tD6 = 1000:100:1e4;
tD_CD = [ tD2 tD3 tD4 tD5 tD6];

% creating various CD and s 
CD = [1  2    4    6  0.2 10  10 .01];
s =  [0 0.5  1.5  2.5  5  4.5 10 60];

for j = 1:length(CD)
    for i = 1:length(tD_CD)
        tD(i) = tD_CD(i)*CD(j);
        % Calculate CDe2s term
        X(j) = CD(j)*exp(2*s(j));
        % Analytical solution of diffusivity equation in Laplace space
        flap = @(p) ( besselk(0, sqrt(p)) + s(j).*sqrt(p).*besselk(1, sqrt(p)) ) ./ ...
            (  p* ( sqrt(p)*besselk(1, sqrt(p)) + CD(j)*p*( besselk(0, sqrt(p)) ...
            + s(j).*sqrt(p).*besselk(1, sqrt(p)) ) )  );
        % Inverting solution numerically using Stehfest algorithm
        PD(i, j) = stehfestAlgorithm(flap, tD(i), 12);
    end
end

% Calculating the derivative plot
derP = zeros(8, length(tD_CD) -1);
for j = 1:length(CD) 
    derP(j, :) = diff(PD(:, j))'.*tD_CD(1:(length(tD_CD)-1))./diff(tD_CD);
end

%% PLOTING THE RESULTS
figure
gg = loglog(tD_CD, PD, 'LineWidth', 2 ); grid on;
set(gg(8),'Color',[0.87058824300766 0.490196079015732 0]);
text(.5,400, 'C_D e^{2s}', 'FontSize',12);

lgd = legend( num2str(round(X(1))) , num2str(round(X(2))), num2str(round(X(3))), ...
    num2str(round(X(4))), num2str(round(X(5))), ['1e', num2str(round(log10(X(6))))], ...
    ['1e', num2str(round(log10((X(7)))) )], ['1e', num2str(round(log10(X(8))))], ...
    'Location', 'NorthWest' );
xlabel('t_D/C_D');
ylabel('Dimensionless Pressure, P_D');
title(' Gringarten Type-Curve');
ylim([1e-1 1e2]);

figure
gg = loglog(tD_CD, PD, 'LineWidth', 2 ); grid on;
set(gg(8),'Color',[0.87058824300766 0.490196079015732 0]);
text(.5,400, 'C_D e^{2s}', 'FontSize',12);

lgd = legend( num2str(round(X(1))) , num2str(round(X(2))), num2str(round(X(3))), ...
    num2str(round(X(4))), num2str(round(X(5))), ['1e', num2str(round(log10(X(6))))], ...
    ['1e', num2str(round(log10((X(7)))) )], ['1e', num2str(round(log10(X(8))))], ...
    'Location', 'NorthWest' );
xlabel('t_D/C_D');
ylabel('Dimensionless Pressure, P_D, t_D/C_D P^\prime_D');
title(' Gringarten Type-Curve');

hold on
gder = loglog(tD_CD(1:(length(tD_CD)-1)), derP', 'LineWidth', 2);
set(gder(8),'Color',[0.87058824300766 0.490196079015732 0]);
grid on;
ylim([1e-1 1e2]);

figure
gder = loglog(tD_CD(1:(length(tD_CD)-1)), derP', 'LineWidth', 2);
set(gder(8),'Color',[0.87058824300766 0.490196079015732 0]);
lgd = legend( num2str(round(X(1))) , num2str(round(X(2))), num2str(round(X(3))), ...
    num2str(round(X(4))), num2str(round(X(5))), ['1e', num2str(round(log10(X(6))))], ...
    ['1e', num2str(round(log10((X(7)))) )], ['1e', num2str(round(log10(X(8))))], ...
    'Location', 'NorthWest' );
xlabel('t_D/C_D');
ylabel('Dimensionless Pressure Derivative, t_D/C_D P^\prime_D');
title(' Bourdet Type-Curve');
grid on;
text(.5,400, 'C_D e^{2s}', 'FontSize',12);
ylim([1e-1 1e2]);