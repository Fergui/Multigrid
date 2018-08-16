function Y=dst2(X);
% Y=dst2(X)
% 2D discrete sine transform
% Jan Mandel
Z=dst(X');
Y=dst(Z');

