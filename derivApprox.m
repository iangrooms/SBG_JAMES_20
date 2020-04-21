function [b] = derivApprox(z, h)
% This function uses a finite volume reconstruction to obtain an approximate derivative of z
% z is a vector of values
% h is a vector of heights

len = length(z) ; 

b = zeros(len, 1) ;

for i = 2:len-1
	denom = ( h(i-1) + h(i) )*( h(i) + h(i+1) )*( h(i-1) + h(i) + h(i+1) ) ;
	num1 = - ( h(i)^2 + 2*h(i+1)^2 + 3*h(i)*h(i+1) ) ; 
	num2 = -2*h(i-1)^2 - 3*h(i-1)*h(i) + 3*h(i)*h(i+1) + 2*h(i+1)^2 ;
	num3 = h(i)^2 + 2*h(i-1)^2 + 3*h(i-1)*h(i) ;
 
	b(i) = ( num1*z(i-1) + num2*z(i) + num3*z(i+1) ) / denom ;
end % end for i

b(1) = b(2) ;
b(len) = b(len-1) ;

end
