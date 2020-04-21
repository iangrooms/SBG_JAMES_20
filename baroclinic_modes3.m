filename = 'current_data';

tic;

ncid = netcdf.open(filename , 'NOWRITE') ;
TAREA = netcdf.getVar(ncid, 17) ; % Area of T-cell
TAREA = TAREA(:,121:1182) ;
UAREA = netcdf.getVar(ncid, 16) ; % Area of U-cell
UAREA = UAREA(:,121:1182) ;
ULAT = netcdf.getVar(ncid, 10) ; % Latitude of U-cell
ULAT = ULAT(:,121:1182) ;
UVEL = netcdf.getVar(ncid, 68) ; % Horizontal velocity in x direction
UVEL = UVEL(:,121:1182,:) ;
VVEL = netcdf.getVar(ncid, 70) ; % Horizontal velocity in y direction
VVEL = VVEL(:,121:1182,:) ;
PD = netcdf.getVar(ncid, 123) ; % Potential Density
PD = PD(:,121:1182,:) ;
KMU = netcdf.getVar(ncid, 14) ; % k index of deepest U grid cell
KMU = KMU(:,121:1182) ;
KMT = netcdf.getVar(ncid,13) ; % k index of deepest T grid cell
KMT = KMT(:,121:1182) ;
DXU = netcdf.getVar(ncid, 20) ;
DXU = DXU(:,121:1182) ;
DYU = netcdf.getVar(ncid, 21) ;
DYU = DYU(:,121:1182) ;
dz = netcdf.getVar(ncid, 7) ; % Thickness of layer k
z_w_bot = netcdf.getVar(ncid, 6) ; % Depth from surface to bottom of layer
netcdf.close(ncid) ;

fprintf('Finished reading netcdf data \n')
toc

KMU = min(KMT,KMU);
clear KMT

f = (2 * 7.292123625e-5).*sind(ULAT) ; % Coriolis frequency
f2 = f.^2 ;

clear( 'ULAT')

% Import smoothed density, U, and V, fields
tic; 
fid = fopen('PDUSmooth3D.dat','r') ;
tmp = fread(fid,1/0,'double') ;
fclose(fid) ; 
tmp=reshape(tmp,[3600 2400 62]) ;
PDSmooth3D = tmp(:,121:1182,:);
clear tmp

fid = fopen('UVELSmooth3D.dat','r') ;
tmp = fread(fid,1/0,'double') ;
fclose(fid) ;
tmp=reshape(tmp,[3600 2400 62]) ;
UVELSmooth3D = tmp(:,121:1182,:);
clear tmp

fid = fopen('VVELSmooth3D.dat','r') ;
tmp = fread(fid,1/0,'double') ;
fclose(fid) ;
tmp=reshape(tmp,[3600 2400 62]) ;
VVELSmooth3D = tmp(:,121:1182,:);
clear tmp
fprintf('Finished reading smoothed data \n')
toc

PDSmooth3D(PDSmooth3D==0) = NaN ;


clear( 'TAREA', 'UAREA' )

PD_INT = zeros(3600, 1062, 61) ; % Calculate PD at interfaces between vertical grid cells
PD_DZ = zeros(3600, 1062, 61) ; % Calculate dPD/dz at interfaces
for d = 1:61
    PD_INT(:, :, d) = ( dz(d+1)/( dz(d) + dz(d+1) ) ) .* PDSmooth3D(:,:, d) + ( dz(d)/( dz(d) + dz(d+1) ) ) .* PDSmooth3D(:, :, d+1)  ;
    PD_DZ(:, :, d) = ( 2 / ( dz(d) + dz(d+1) ) ) .* ( PDSmooth3D(:, :, d) - PDSmooth3D(:, :, d+1) ) ;
end

N2 = (-980.6 .* PD_DZ) ./ PD_INT ; % Buoyancy frequency
N2(N2 < 1e-10) = 1e-10 ; 
N = sqrt(N2);

S = bsxfun(@times,f2, 1./ N2) ; % S(z) = f^2 / N^2(z)

PD_DX = zeros(3600, 1062, 62) ;
PD_DY = zeros(3600, 1062, 62) ;

for z = 1:62
    for y = 1:1062 
        PD_DX(:,y,z) = derivApprox( PDSmooth3D(:,y,z), DXU(:,y) ) ;
    end
    for x = 1:3600
        PD_DY(x, :, z) = derivApprox( PDSmooth3D(x,:,z), DYU(x,:) ) ;
    end
end 

PD_GRAD = sqrt(PD_DX.^2 + PD_DY.^2) ; % This has spurious values next to land
PD_GRAD(isnan(PD_GRAD)) = 0;

fprintf('Finished computing derivatives of mean density \n')
toc

clear( 'DXU', 'DYU' , 'PD_DX', 'PD_DY')

LD = zeros(3600, 1062) ; % length scale (deformation radius)
UZ = zeros(3600, 1062) ; % projection of v onto z
VZ = zeros(3600, 1062) ; % projection of u onto z
NORM_UV = zeros(3600, 1062) ; % norm(u)^2 + norm(v)^2
PD_COEF1 = zeros(3600, 1062) ; % EBT projection coeficients
PD_COEF2 = zeros(3600, 1062) ; % EBT x grad projection coeficients
NORM_PD = zeros(3600,1062); % norm (pd)^2

% Calculate at all locations: relative residual norm, projection of u onto z
% projection of v onto z, deformation radius (related to first eigenvalue),
% dimensional time scale

for y = 1:1062
for x = 1:3600

% Calculate first eigenvector (z)
depth = KMU(x, y) - 2 ;

if depth > 15

% Construct tri-diagonal matrix which is the discretization of the PV stretching operator
b = zeros(depth, 1) ;
a = zeros(depth - 1, 1) ;
c = zeros(depth - 1, 1) ;

SLoc = S(x,y,:);
SLoc = SLoc(:);

b(1) = ( -2 * SLoc(1) ) / ( dz(1) * ( dz(1)+dz(2) ) ) ;
c(1) = -b(1);
for i = 2:depth - 1
        c(i) = ( 2 * SLoc(i) ) / ( dz(i) * ( dz(i) + dz(i+1) ) ) ;
        a(i-1) = ( 2 * SLoc(i-1) ) / ( dz(i) * ( dz(i) + dz(i-1) ) ) ;
        b(i) =  -a(i-1)-c(i) ;
end
a(depth-1) = 2*SLoc(depth-1)/( dz(depth)*(dz(depth)+dz(depth-1)) );
b(depth) = (-2*SLoc(depth-1)/dz(depth))*(1/(dz(depth)+dz(depth-1)) + 1/dz(depth) );

% Apply 5 steps of the inverse power method to approximate the first eigenvector
z = double(0:(depth-1)) ;
for i = 1:5
    z = tridisolve(a, b, c, z) ;
    z = z / sqrt(dot(z, z)) ;
end

% Compute eddy length scale
% Calculate deformation radius as -1/l1 = zT A z  (Rayleigh quotient)
on_diag = dot(b, z.^2) ;
z_diff = z(1:depth-1) .* z(2:depth) ;
off_diag = dot(a+c, z_diff) ;
quad_form = on_diag + off_diag ;
LD(x,y) = sqrt( -1 / quad_form ) /1E5; % Deformation Radius in km

% Normalize the eigenvector
z_flip = z;%flip(z) ;
normz = baro_inner_product(z_flip, z_flip, dz) ;
z_flip = z_flip ./ sqrt(normz) ;

% Extract eddy velocity
u = UVEL(x,y, 1:depth) - UVELSmooth3D(x, y, 1:depth) ;
v = VVEL(x,y, 1:depth) - VVELSmooth3D(x, y, 1:depth) ;
u = u(:) ;
v = v(:) ;

% Project eddy velocity onto eigenvector
uDotz = baro_inner_product(u, z_flip, dz) ;
vDotz = baro_inner_product(v, z_flip, dz) ;
UZ(x, y) = uDotz ;     % projection of u onto z
VZ(x, y) = vDotz ;     % projection of v onto z

% Compute EKE
normu2 = baro_inner_product(u, u, dz) ;
normv2 = baro_inner_product(v, v, dz) ;
NORM_UV(x, y) = normu2 + normv2 ; % 2*(depth-averaged EKE)

% Extract eddy density
pd = PD(x, y, 1:depth) - PDSmooth3D(x, y, 1:depth) ;
pd = pd(:) ;

% (1) Vertical density on EBC Mode
PD_COEF1(x,y) = baro_inner_product(pd, z_flip, dz);
NORM_PD(x,y) = baro_inner_product(pd, pd, dz);

% (2) Vertical density structure proportional to EBC Mode * grad rho
pd_grad = PD_GRAD(x,y,1:depth) ;
pd_grad = pd_grad(:) ;
z_flip = z_flip(:) ;
z_vert1 = pd_grad(1:depth) .* z_flip ;
normz_vert1 = baro_inner_product(z_vert1, z_vert1, dz) ;
z_vert1 = z_vert1 ./ sqrt(normz_vert1) ;
PD_COEF2(x,y) = baro_inner_product(pd,z_vert1,dz);
end % end if

end % end for x
end % end for y

toc

fid = fopen('LD.dat','w');fwrite(fid,LD,'double');fclose(fid);
fid = fopen('UZ.dat','w');fwrite(fid,UZ,'double');fclose(fid);
fid = fopen('VZ.dat','w');fwrite(fid,VZ,'double');fclose(fid);
fid = fopen('NORM_UV.dat','w');fwrite(fid,NORM_UV,'double');fclose(fid);
fid = fopen('PD_COEF1.dat','w');fwrite(fid,PD_COEF1,'double');fclose(fid);
fid = fopen('PD_COEF2.dat','w');fwrite(fid,PD_COEF2,'double');fclose(fid);
fid = fopen('NORM_PD.dat','w');fwrite(fid,NORM_PD,'double');fclose(fid);

toc
