function[prod] = baro_inner_product(a, b, dz)

depth = max(size(a)) ;
prod = 0 ;
sumdz = 0 ;

for i=1:depth
	prod = prod + (a(i)*b(i)*dz(i)) ;
	sumdz = sumdz + dz(i) ;
end

prod = prod / sumdz ;

end
