filename = 'current_data'; % A netcdf data file containing a 5-day-average output from the Johnson et al. (2016) simulation

%% Read netcdf data
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

KMU = min(KMT,KMU); % Depth based on the smaller of T or U points
clear KMT

f = (2 * 7.292123625e-5).*sind(ULAT) ; % Coriolis frequency

clear( 'ULAT')

% Import smoothed fields, output from smooth.f90
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
LD = zeros(3600, 1062) ; % length scale (deformation radius)

NORM_UV = zeros(3600, 1062) ; % norm(u)^2 + norm(v)^2
UZ = zeros(3600, 1062) ; % projection of v onto z
VZ = zeros(3600, 1062) ; % projection of u onto z

NORM_PD = zeros(3600,1062); % variance of PD
APE = zeros(3600, 1062); % (b'/N)^2 averaged over depth
PD_COEF3 = zeros(3600, 1062) ; % EBT*N^2 projection coefficients
PD_COEF2 = zeros(3600, 1062, 2) ; % EBT*GRAD projection coefficients
PD_COEF4 = zeros(3600, 1062, 2) ; % EBT*EOF projection coefficients
PVE = zeros(3600,1062); % Percent Variance Explained by leading EOF

NORM_FR = zeros(3600,1062); % Norm of SGS density flux
FX_COEF1 = zeros(3600,1062,2); % EBT^2*GRAD projection coefficients
FY_COEF1 = zeros(3600,1062,2); % EBT^2*GRAD projection coefficients
FX_COEF2 = zeros(3600,1062); % EBT*GRAD projection coefficients
FY_COEF2 = zeros(3600,1062); % EBT*GRAD projection coefficients
FX_COEF3 = zeros(3600,1062); % N^2*GRAD projection coefficients
FY_COEF3 = zeros(3600,1062); % N^2*GRAD projection coefficients


%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

for y = 1:1062
for x = 1:3600

% Calculate first eigenvector (z)
depth = KMU(x, y) - 3 ;

if depth > 42 % 2 km total depth, but not using bottom 3

% Te
M2 = 980.6*squeeze(PD_GRAD(x,y,1:depth)./PDBar(x,y,1:depth));
Te(x,y) = sum(sqrt(squeeze(N2_v(x,y,1:depth))).*dz(1:depth))/sum(M2.*dz(1:depth)); % [s]

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

% Compute eddy length scale
% Calculate deformation radius as -1/l1 = zT A z  (Rayleigh quotient)
on_diag = dot(b, z.^2) ;
z_diff = z(1:depth-1) .* z(2:depth) ;
off_diag = dot(a+c, z_diff) ;
quad_form = on_diag + off_diag ;
LD(x,y) = sqrt( -1 / quad_form ) /1E5; % Deformation Radius [km]

% Normalize eigenvector, it is not nondimensional
z = z / sqrt(baro_inner_product(z,z,dz));

% Extract eddy velocity
u = squeeze(UVEL(x,y, 1:depth));
v = squeeze(VVEL(x,y, 1:depth));

% Project eddy velocity onto eigenvector
UZ(x, y) = baro_inner_product(u, z, dz); % projection of u onto z [cm/s]
VZ(x, y) = baro_inner_product(v, z, dz); % projection of v onto z [cm/s]

% Compute EKE
normu2 = baro_inner_product(u, u, dz) ;
normv2 = baro_inner_product(v, v, dz) ;
NORM_UV(x, y) = normu2 + normv2 ; % 2*(depth-averaged EKE) [cm^2/s^2]

% Extract eddy density
pd = squeeze(PD(x, y, 1:depth) - PDBar(x, y, 1:depth));

% Eddy density variance
NORM_PD(x,y) = baro_inner_product(pd,pd,dz); 

% Eddy APE
pd_b = -980.6*pd./squeeze(PDBar(x,y,1:depth));
APE(x,y) = baro_inner_product(pd_b,pd_b./max(squeeze(N2_v(x,y,1:depth)),1E-8),dz);

% (1) Vertical density on EBC Mode * N^2
z_vert1 = squeeze(N2_v(x,y,1:depth)).*z ;
z_vert1 = z_vert1 / sqrt(baro_inner_product(z_vert1, z_vert1, dz)) ; % normalize
PD_COEF3(x,y) = baro_inner_product(pd, z_vert1, dz); % [g/cm^3]

% (2) Vertical density structure proportional to EBC Mode * |grad rho|
% Normalized
pd_grad = squeeze(PD_GRAD(x,y,1:depth));
z_vert1 = pd_grad.*z ; % Not normalized
tmp = baro_inner_product(z_vert1,z_vert1,dz);
PD_COEF2(x,y,1) = baro_inner_product(pd,z_vert1,dz) / tmp; % [cm]
PD_COEF2(x,y,2) = PD_COEF2(x,y,1)*sqrt(tmp); % [g/cm^3]

% (3) Construct EBT*PD_DX and EBT*PD_DY, scale with sqrt(dz/H)
A = diag(sqrt(dz(1:depth)/sum(dz(1:depth))).*z)*[squeeze(PD_DX(x,y,1:depth)) squeeze(PD_DY(x,y,1:depth))];
if(~any(isnan(A(:))))
[P,Sigma,Q] = svd(A);
PVE(x,y) = Sigma(1,1)^2/(Sigma(1,1)^2+Sigma(2,2)^2);
z_vert1 = sqrt(sum(dz(1:depth)))*Sigma(1,1)*P(:,1)./sqrt(dz(1:depth));
z_vert1 = -z_vert1*sign(z_vert1(1)); % Sign convention on EOF

% Project onto leading EOF
tmp = baro_inner_product(z_vert1,z_vert1,dz);
PD_COEF4(x,y,1) = baro_inner_product(pd,z_vert1,dz) / tmp; % [cm]
PD_COEF4(x,y,2) = PD_COEF4(x,y,1)*sqrt(tmp); % [g/cm^3]
end % end if nan

% (4) SGS density flux: EBT^2*GRAD
fx = squeeze(FX(x,y,1:depth));
fy = squeeze(FY(x,y,1:depth));
normu2 = baro_inner_product(fx, fx, dz) ;
normv2 = baro_inner_product(fy, fy, dz) ;
NORM_FR(x, y) = normu2 + normv2 ; % depth-average norm of SGS density flux

z_vert1 = pd_grad(1:depth).*(z.^2); 
tmp = baro_inner_product(z_vert1, z_vert1, dz);
FX_COEF1(x,y,1) = baro_inner_product(fx, z_vert1, dz) / tmp; % [cm^2/s]
FY_COEF1(x,y,1) = baro_inner_product(fy, z_vert1, dz) / tmp; % [cm^2/s]

FX_COEF1(x,y,2) = FX_COEF1(x,y,1)*sqrt(tmp); % [g / cm^2 s]
FY_COEF1(x,y,2) = FY_COEF1(x,y,1)*sqrt(tmp); % [g / cm^2 s]

% (5) SGS density flux: EBT*GRAD
z_vert1 = pd_grad(1:depth).*z; 
z_vert1 = z_vert1 / sqrt(baro_inner_product(z_vert1, z_vert1, dz)) ; % normalize
FX_COEF2(x,y) = baro_inner_product(fx, z_vert1, dz); % [cm^2/s]
FY_COEF2(x,y) = baro_inner_product(fy, z_vert1, dz); % [cm^2/s]

% (6) SGS density flux: N^2*GRAD
z_vert1 = pd_grad(1:depth).*squeeze(N2_v(x,y,1:depth)); 
z_vert1 = z_vert1 / sqrt(baro_inner_product(z_vert1, z_vert1, dz)) ; % normalize
FX_COEF3(x,y) = baro_inner_product(fx, z_vert1, dz); % [cm^2/s]
FY_COEF3(x,y) = baro_inner_product(fy, z_vert1, dz); % [cm^2/s]

end % end if

end % end for x
end % end for y
fprintf('Finished loop over x & y \n')
toc

save('VertStruct.mat','LD','Te','NORM_UV','UZ','VZ',...
'NORM_PD','PD_COEF2','PD_COEF3','PD_COEF4','APE','NORM_FR','FX_COEF1','FY_COEF1',...
'FX_COEF2','FY_COEF2','FX_COEF3','FY_COEF3','-v7.3')
