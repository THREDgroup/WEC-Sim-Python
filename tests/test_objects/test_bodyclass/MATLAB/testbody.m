clc; clear all; close all;
A = zeros(2,1,2);
%A(:,:,1) = [1 2]';
A(:,:,2) = [3 4]';
A
B = squeeze(A)