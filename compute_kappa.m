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
N2_v(:,:,62) = NaN;
N2_i(:,:,61) = NaN;
N2_v = max(N2_v,1e-10);
N2_i = max(N2_i,1e-10);
S = bsxfun(@times,f.^2, 1./ N2_i) ; % S(z) = f^2 / N^2(z)
clear N2_i f
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

fprintf(' Horizontal Derivatives \n')
toc

Err_K = zeros(3600,1062,3); % EKE error
K_X = zeros(3600,1062,6); 
K_Y = zeros(3600,1062,6); 
K_C = zeros(3600,1062,6); 
Err_F = zeros(3600,1062,6); % square error in projection

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

for y = 1:1062
for x = 1:3600

% Calculate first eigenvector (z)
depth = KMU(x, y) - 3 ;

if depth > 42 % 2 km total depth, but not using bottom 3

% Construct tri-diagonal matrix which is the discretization of the PV
% stretching operator
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

% Normalize eigenvector, it is not nondimensional yet
z = z / sqrt(baro_inner_product(z,z,dz));

%% Eddy velocity

% Extract eddy velocity
u = squeeze(UVEL(x,y, 1:depth));
v = squeeze(VVEL(x,y, 1:depth));

w = diag(dz(1:depth)/sum(dz(1:depth)));
BT = ones(depth,1); % BT mode

% (1) Project onto EBT
uErr = u - z*(z'*w*u);
vErr = v - z*(z'*w*v);
Err_K(x,y,1) = uErr'*w*uErr + vErr'*w*vErr;

% (2) Project onto BT
uErr = u - BT*(BT'*w*u);
vErr = v - BT*(BT'*w*v);
Err_K(x,y,2) = uErr'*w*uErr + vErr'*w*vErr;

% (3) EBT & BT
w = diag(sqrt(dz(1:depth)/sum(dz(1:depth))));
Vec = [z BT];
uTmp = (w*Vec)\(w*u);
vTmp = (w*Vec)\(w*v);
uErr = w*(u - Vec*uTmp);
vErr = w*(v - Vec*vTmp);
Err_K(x,y,3) = uErr'*uErr + vErr'*vErr;

%% SGS density flux
fx = squeeze(FX(x,y,1:depth));
fy = squeeze(FY(x,y,1:depth));

% Density gradient
pdx = squeeze(PD_DX(x,y,1:depth));
pdy = squeeze(PD_DY(x,y,1:depth));

% (1) Projection using iso kappa, depth-indep
w = diag(dz(1:depth)/sum(dz(1:depth)));
K_X(x,y,1) = -(fx'*w*pdx + fy'*w*pdy)/(pdx'*w*pdx + pdy'*w*pdy);
K_Y(x,y,1) = K_X(x,y,1);
K_C(x,y,1) = 0;

Err_F(x,y,1) = (fx + K_X(x,y,1)*pdx)'*w*(fx + K_X(x,y,1)*pdx) + ...
    (fy + K_X(x,y,1)*pdy)'*w*(fy + K_X(x,y,1)*pdy);

% (2) Projection using aniso kappa, depth-indep
w = diag(sqrt([dz(1:depth);dz(1:depth)]/sum(dz(1:depth))));
A = zeros(2*depth,3);
A(1:depth,1) = -pdx;
A(1:depth,3) = -pdy;
A(depth+1:2*depth,2) = -pdy;
A(depth+1:2*depth,3) = -pdx;
A = w*A;
b = w*[fx;fy];

tmp = A\b;
K_X(x,y,2) = tmp(1);
K_Y(x,y,2) = tmp(2);
K_C(x,y,2) = tmp(3);

Err_F(x,y,2) = norm(b-A*tmp,2)^2;

% (3) Projection using iso kappa ~ phi
w = diag(dz(1:depth)/sum(dz(1:depth)));
zx = z.*pdx;
zy = z.*pdy;
K_X(x,y,3) = -(fx'*w*zx + fy'*w*zy)/(zx'*w*zx + zy'*w*zy);
K_Y(x,y,3) = K_X(x,y,3);
K_C(x,y,3) = 0;

Err_F(x,y,3) = (fx + K_X(x,y,3)*zx)'*w*(fx + K_X(x,y,3)*zx) + ...
    (fy + K_X(x,y,3)*zy)'*w*(fy + K_X(x,y,3)*zy);

% (4) Projection using aniso kappa ~ phi
w = diag(sqrt([dz(1:depth);dz(1:depth)]/sum(dz(1:depth))));
A = zeros(2*depth,3);
A(1:depth,1) = -zx;
A(1:depth,3) = -zy;
A(depth+1:2*depth,2) = -zy;
A(depth+1:2*depth,3) = -zx;
A = w*A;
b = w*[fx;fy];

tmp = A\b;
K_X(x,y,4) = tmp(1);
K_Y(x,y,4) = tmp(2);
K_C(x,y,4) = tmp(3);

Err_F(x,y,4) = norm(b-A*tmp,2)^2;

% (5) Projection using iso kappa ~ N^2
w = diag(dz(1:depth)/sum(dz(1:depth)));
N2Prof = squeeze(N2_v(x,y,1:depth));
N2Prof = max(N2Prof/max(N2Prof),0.1);
zx = pdx.*N2Prof;
zy = pdy.*N2Prof;
K_X(x,y,5) = -(fx'*w*zx + fy'*w*zy)/(zx'*w*zx + zy'*w*zy);
K_Y(x,y,5) = K_X(x,y,5);
K_C(x,y,5) = 0;

Err_F(x,y,5) = (fx + K_X(x,y,5)*zx)'*w*(fx + K_X(x,y,5)*zx) + ...
    (fy + K_X(x,y,5)*zy)'*w*(fy + K_X(x,y,5)*zy);

% (6) Projection using aniso kappa~N^2
w = diag(sqrt([dz(1:depth);dz(1:depth)]/sum(dz(1:depth))));
A = zeros(2*depth,3);
A(1:depth,1) = -zx;
A(1:depth,3) = -zy;
A(depth+1:2*depth,2) = -zy;
A(depth+1:2*depth,3) = -zx;
A = w*A;
b = w*[fx;fy];

tmp = A\b;
K_X(x,y,6) = tmp(1);
K_Y(x,y,6) = tmp(2);
K_C(x,y,6) = tmp(3);

Err_F(x,y,6) = norm(b-A*tmp,2)^2;

end % end if

end % end for x
end % end for y
fprintf('Finished loop over x & y \n')
toc

save('Kappa_EKE.mat','K_X','K_Y','K_C','Err_F','Err_K','-v7.3')
