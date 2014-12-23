function [] = Plot_EV_CB(L_A,L_J,Slopes,Alpha)

if nargin<1
    Data = load('Gen4CLSaltMat.mat');
    L_A = Data.AnEigV;
    L_J = Data.NumEigV;
    Slopes = Data.Slopes;
    Alpha = Data.Alpha;
end

%%%%%%%% Lambda vs Slope Comparison %%%%%%%
figure()
h = plot(Slopes(1:1:end),abs(L_A(1,1:1:end))','-*b');hold on;
plot(Slopes(1:1:end),abs(L_A)','-*')
h1 = plot(Slopes(1:1:end),abs(L_J(1,1:1:end)'),'-ob');
plot(Slopes(1:1:end),abs(L_J(1:5,:))','-o');
legend([h,h1],{'Numeric LC + Analytic EV','Numeric (LC+EV)'}),%'Numeric LC+Analytic+PWL EV'
xlabel( 'Slope \gamma^\circ', 'FontSize',15,'FontWeight','bold' )
ylabel( '\lambda', 'FontSize',15,'FontWeight','bold')

%%%%%%%% Lambda vs Alpha Comparison %%%%%%%
figure()
h = plot(Alpha(1:1:end),abs(L_A(1,1:1:end))','-*b');hold on;
plot(Alpha(1:1:end),abs(L_A)','-*')
h1 = plot(Alpha(1:1:end),abs(L_J(1,1:1:end)'),'-ob');
plot(Alpha(1:1:end),abs(L_J(1:5,:))','-o');

% plot(Alpha(1:1:end),abs(L_A(1,1:1:end))','-*b')
% hold on;
% h = plot(Alpha(1:1:end),abs(L_A(1,1:1:end))','-*b');
% plot(Alpha(1:1:end),abs(L_A(2,1:1:end))','-*r')
% plot(Alpha(1:1:end),abs(L_A(3,1:1:end))','-*g')%m
% plot(Alpha(1:1:end),abs(L_A(4,1:1:end))','-*m')%k
% plot(Alpha(1:1:end),abs(L_A(5,1:1:end))','-*k')%g
% 
% h1 = plot(Alpha(1:1:end),abs(L_J(1,1:1:end)'),'-ob');
% plot(Alpha(1:1:end),abs(L_J(1,1:1:end)'),'-ob')
% plot(Alpha(1:1:end),abs(L_J(2,1:1:end)'),'-or')
% plot(Alpha(1:1:end),abs(L_J(3,1:1:end)'),'-og')
% plot(Alpha(1:1:end),abs(L_J(4,1:1:end)'),'-om')
% plot(Alpha(1:1:end),abs(L_J(5,1:1:end)'),'-ok')
legend([h,h1],{'Numeric LC + Analytic EV','Numeric (LC+EV)'}),%'Numeric LC+Analytic+PWL EV'
xlabel( '\alpha [rad]', 'FontSize',15,'FontWeight','bold' )
ylabel( '\lambda', 'FontSize',15,'FontWeight','bold')

%%%%%%%% R Locus Comparison %%%%%%%
figure()
Circle();
hold on
plot(real( L_A(1,:)),imag(L_A(1,:)),'.b')%k
plot(real( L_A(2,:)),imag(L_A(2,:)),'.r')%b
plot(real( L_A(3,:)),imag(L_A(3,:)),'.k')%r
plot(real( L_A(4,:)),imag(L_A(4,:)),'.m')%m
plot(real( L_A(5,:)),imag(L_A(5,:)),'.g')%g

plot(real( L_J(1,:)),imag(L_J(1,:)),'ob')
plot(real( L_J(2,:)),imag(L_J(2,:)),'or')
plot(real( L_J(3,:)),imag(L_J(3,:)),'ok')
plot(real( L_J(4,:)),imag(L_J(4,:)),'om')
plot(real( L_J(5,:)),imag(L_J(5,:)),'og')
title('Eigenvalue comparisson')

%%%%%%%% Alpha vs Slope %%%%%%%
figure()
plot(Slopes,Alpha)
hold on
ylabel( '\alpha [rad]', 'FontSize',15,'FontWeight','bold' )
xlabel( 'Slope \gamma^\circ ', 'FontSize',15,'FontWeight','bold' )