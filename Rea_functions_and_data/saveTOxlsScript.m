
close all; clear all; clc

load('MatsRandomRes_4_11_2016.mat','nSims','results','periods')

check_period =  isnan(periods(1,:));
N = sum(~check_period);
v = find(~check_period);

dataMatrix = zeros(19,N);

for i=1:N
    
    dataMatrix(i,1) = results(v(i).seq(1));         % tau
    dataMatrix(i,2) = results(v(i)).seq(2);           % b
    dataMatrix(i,3:6) = (results(v(i)).seq(3:6))';      % c
    c = (results(v(i)).seq(3:6))';
    a = results(v(i)).seq(7:18);
    Aorig = [0    ,a(1) ,a(2) ,a(3);
        a(4) ,0    ,a(5) ,a(6);
        a(7) ,a(8) ,0    ,a(9);
        a(10),a(11),a(12),0   ];
    A = (diag(1./c)*(diag(c)*Aorig)')';
    
    w12 = A(1,2);     w13 = A(1,3);     w14 = A(1,4);
    w21 = A(2,1);     w23 = A(2,3);     w24 = A(2,4);
    w31 = A(3,1);     w32 = A(3,2);     w34 = A(3,4);
    w41 = A(4,1);     w42 = A(4,2);     w43 = A(4,3);
    
    dataMatrix(i,7) = w12;   dataMatrix(i,8) = w12;   dataMatrix(i,9) = w12;
    dataMatrix(i,10) = w21;   dataMatrix(i,11) = w23;   dataMatrix(i,12) = w24;
    dataMatrix(i,13) = w31;   dataMatrix(i,14) = w32;   dataMatrix(i,15) = w34;
    dataMatrix(i,16) = w41;   dataMatrix(i,17) = w42;   dataMatrix(i,18) = w43;

    dataMatrix(i,19) = periods(v(i),1);
end

filename = 'dataMatrixXLSX.xlsx';
xlswrite(filename,dataMatrix);