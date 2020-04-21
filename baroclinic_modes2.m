filename = 'current_data';

tic;

ncid = netcdf.open(filename , 'NOWRITE') ;
ULAT = double(netcdf.getVar(ncid, 10)) ; % Latitude of U-cell
ULAT = ULAT(:,121:1182) ;
UVEL = double(netcdf.getVar(ncid, 68)) ; % Horizontal velocity in x direction
UVEL = UVEL(:,121:1182,:) ;
VVEL = double(netcdf.getVar(ncid, 70)) ; % Horizontal velocity in y direction
VVEL = VVEL(:,121:1182,:) ;
PD = double(netcdf.getVar(ncid, 123)) ; % Potential Density
PD = PD(:,121:1182,:) ;
KMU = netcdf.getVar(ncid, 14) ; % k index of deepest U grid cell
KMU = KMU(:,121:1182) ;
KMT = netcdf.getVar(ncid,13) ; % k index of deepest T grid cell
KMT = KMT(:,121:1182) ;
DXU = double(netcdf.getVar(ncid, 20)) ;
DXU = DXU(:,121:1182) ;
DYU = double(netcdf.getVar(ncid, 21)) ;
DYU = DYU(:,121:1182) ;
dz = double(netcdf.getVar(ncid, 7)) ; % Thickness of layer k
netcdf.close(ncid) ;

fprintf('Finished reading netcdf data \n')
toc

KMU = min(KMT,KMU);
clear KMT

f = (2 * 7.292123625e-5).*sind(ULAT) ; % Coriolis frequency

clear( 'ULAT')

% Import smoothed density field
tic;
fid = fopen('UVELSmooth3D.dat','r') ;
tmp = fread(fid,3600*2400*62,'double') ;
fclose(fid) ;

tmp = reshape(tmp,[3600 2400 62]);
UBar = tmp(:,121:1182,:);
clear tmp
UVEL = UVEL - UBar; % Eddy velocity

fid = fopen('VVELSmooth3D.dat','r') ;
tmp = fread(fid,3600*2400*62,'double') ;
fclose(fid) ;

tmp = reshape(tmp,[3600 2400 62]);
VBar = tmp(:,121:1182,:);
clear tmp
VVEL = VVEL - VBar; % Eddy velocity

fid = fopen('PDUSmooth3D.dat','r') ;
tmp = fread(fid,3600*2400*62,'double') ;
fclose(fid) ;

tmp = reshape(tmp,[3600 2400 62]);
PDBar = tmp(:,121:1182,:);
clear tmp

fid = fopen('UPDSmooth3D.dat','r') ;
tmp = fread(fid,3600*2400*62,'double') ;
fclose(fid) ;

tmp = reshape(tmp,[3600 2400 62]);
FX = tmp(:,121:1182,:);
clear tmp
FX = FX - UBar.*PDBar;
clear UBar

fid = fopen('VPDSmooth3D.dat','r') ;
tmp = fread(fid,3600*2400*62,'double') ;
fclose(fid) ;

tmp = reshape(tmp,[3600 2400 62]);
FY = tmp(:,121:1182,:);
clear tmp
FY = FY - VBar.*PDBar;
clear VBar

fprintf('Finished reading smoothed data \n')
toc

PDBar(PDBar==0) = NaN ;

N2_v = zeros(3600, 1062, 62) ; % -(g/rho)*drho/dz
N2_i = zeros(3600,1062,61) ; % value at interface
N2_v(:,:,1) = dz(1)^2*(-PDBar(:,:,2)+PDBar(:,:,3)) + 2*dz(3)^2*(PDBar(:,:,1)-PDBar(:,:,2)) ...
                + 2*dz(2)^2*(2*PDBar(:,:,1)+PDBar(:,:,3)-3*PDBar(:,:,2));
N2_v(:,:,1) = N2_v(:,:,1) + (6*dz(2)*dz(3) + 3*dz(1)*dz(3))*(PDBar(:,:,1)-PDBar(:,:,2)) ...
              +3*dz(2)*dz(1)*(PDBar(:,:,1)-2*PDBar(:,:,2)+PDBar(:,:,3));
N2_v(:,:,1) = N2_v(:,:,1) ./ ((dz(2)+dz(1))*(dz(2)+dz(3))*(dz(3)+dz(2)+dz(1)));
N2_v(:,:,1) = -980.6*N2_v(:,:,1)./PDBar(:,:,1);
for d = 2:61
% Code below is based on piecewise-parabolic FV reconstruction, returns average of N^2 over cell
    N2_v(:,:,d) = 2*dz(d-1)^2*(PDBar(:,:,d)-PDBar(:,:,d+1)) + 2*dz(d+1)^2*(PDBar(:,:,d-1)-PDBar(:,:,d)) ...
                + dz(d)^2*(PDBar(:,:,d-1)-PDBar(:,:,d+1)) ...
              + 3*dz(d)*(dz(d-1)*(PDBar(:,:,d)-PDBar(:,:,d+1)) + dz(d+1)*(PDBar(:,:,d-1)-PDBar(d)));
    N2_v(:,:,d) = N2_v(:,:,d) / ((dz(d)+dz(d-1))*(dz(d)+dz(d+1))*(dz(d+1)+dz(d)+dz(d-1)));
    N2_v(:,:,d) = -980.6*N2_v(:,:,d) ./ PDBar(:,:,d);
% Code below is based on piecewise-linear FV reconstruction, returns N^2 at interfaces. Very noisy & similar to above.
    N2_i(:,:,d-1) = -980.6*2*(PDBar(:,:,d-1)-PDBar(:,:,d))./(dz(d)*PDBar(:,:,d-1) + dz(d-1)*PDBar(:,:,d));
end
N2_v(:,:,62) = NaN;
N2_i(:,:,61) = NaN;
N2_v = max(N2_v,1e-10);
N2_i = max(N2_i,1e-10);
S = bsxfun(@times,f.^2, 1./ N2_i) ; % S(z) = f^2 / N^2(z)
clear N2_i
fprintf(' Finished computing N^2 and S \n')
toc

for z = 1:62
    for y = 1:1062
        PD_DX(:,y,z) = derivApprox( PDBar(:,y,z), DXU(:,y) ) ;
    end
    for x = 1:3600
        PD_DY(x, :, z) = derivApprox( PDBar(x,:,z), DYU(x,:) ) ;
    end
end
clear DXU DYU

% Set spurious values next to land to zero
PD_DX(isnan(PD_DX)) = 0;
PD_DY(isnan(PD_DY)) = 0;

PD_GRAD = sqrt(PD_DX.^2 + PD_DY.^2) ;
fprintf(' Horizontal Derivatives \n')
toc

Te = zeros(3600,1062); % N/M^2

U0 = zeros(3600, 1062) ; % projection of v onto bt mode
V0 = zeros(3600, 1062) ; % projection of u onto bt mode

APE = zeros(3600, 1062); % (b'/N)^2 averaged over depth
PD_COEF5 = zeros(3600, 1062, 2) ; % GRAD projection coefficients

FX_COEF2 = zeros(3600,1062,2); % EBT*GRAD projection coefficients
FY_COEF2 = zeros(3600,1062,2); % EBT*GRAD projection coefficients


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

for y = 1:1062
for x = 1:3600

% Calculate first eigenvector (z)
depth = KMU(x, y) - 3 ;

if depth > 42 % 2 km total depth, but not using bottom 3

% Te
M2 = 980.6*squeeze(PD_GRAD(x,y,1:depth)./PDBar(x,y,1:depth));
Te(x,y) = sum(sqrt(max(squeeze(N2_v(x,y,1:depth)),1E-8)).*dz(1:depth))/sum(M2.*dz(1:depth)); % [s]

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

% Apply 15 steps of the inverse power method to approximate the first eigenvector
z = double(0:(depth-1))' ;
for i = 1:15
    z = tridisolve(a, b, c, z) ;
    z = z / sqrt(z'*z) ;
end
z = z(:)*sign(z(1));

% Normalize eigenvector, it is not nondimensional
z = z / sqrt(baro_inner_product(z,z,dz));

% Extract eddy velocity
u = squeeze(UVEL(x,y, 1:depth));
v = squeeze(VVEL(x,y, 1:depth));

% Project eddy velocity onto barotropic mode 
U0(x, y) = baro_inner_product(u, ones(depth,1), dz); % projection of u onto z [cm/s]
V0(x, y) = baro_inner_product(v, ones(depth,1), dz); % projection of v onto z [cm/s]

% Extract eddy density
pd = squeeze(PD(x, y, 1:depth) - PDBar(x, y, 1:depth));

% Eddy APE
pd_b = -980.6*pd./squeeze(PDBar(x,y,1:depth));
APE(x,y) = baro_inner_product(pd_b,pd_b./max(squeeze(N2_v(x,y,1:depth)),1E-6),dz);

% (*) Vertical density structure proportional to |grad rho|
pd_grad = squeeze(PD_GRAD(x,y,1:depth));
z_vert1 = pd_grad ; % Not normalized
tmp = baro_inner_product(z_vert1,z_vert1,dz);
PD_COEF5(x,y,1) = baro_inner_product(pd,z_vert1,dz) / tmp; % [cm]
PD_COEF5(x,y,2) = PD_COEF5(x,y,1)*sqrt(tmp); % [g/cm^3]

% SGS density flux
fx = squeeze(FX(x,y,1:depth));
fy = squeeze(FY(x,y,1:depth));

% (5) SGS density flux: EBT*GRAD
z_vert1 = pd_grad(1:depth).*z; 
tmp = baro_inner_product(z_vert1, z_vert1, dz);
FX_COEF2(x,y,1) = baro_inner_product(fx, z_vert1, dz) / tmp; % [cm^2/s]
FY_COEF2(x,y,1) = baro_inner_product(fy, z_vert1, dz) / tmp; % [cm^2/s]

FX_COEF2(x,y,2) = FX_COEF2(x,y,1)*sqrt(tmp); % [g / cm^2 s]
FY_COEF2(x,y,2) = FY_COEF2(x,y,1)*sqrt(tmp); % [g / cm^2 s]

end % end if

end % end for x
end % end for y
fprintf('Finished loop over x & y \n')
toc

save('VertStruct2.mat','Te','APE','U0','V0','PD_COEF5','FX_COEF2','FY_COEF2','-v7.3')
