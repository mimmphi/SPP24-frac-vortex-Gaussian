function [Q,T,F] = linfoots(I_o,I_r)
%I_o: original image
%I_r: reconstructed image
Q = sum(sum(I_o.*I_r))/sum(sum(I_o.*I_o));
T = sum(sum(I_r.*I_r))/sum(sum(I_o.*I_o));
F = 2*Q - T;
end