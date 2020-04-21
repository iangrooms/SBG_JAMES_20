load Depth_Land_Masks.mat mask % the mask array has 1 in locations with > 2km depth and >200 points from shore
cd 1Deg/
d = dir;
cd ../
TOL = 2.5E-9; % (g/cm^3)^2
NPD = [];
PD2 = []; % EBT*GRAD
%PD3 = []; % EBT*N
%PD4 = []; % EBT*EOF
PD5 = []; % GRAD
ER6 = []; % {phi,psi_0}GRAD
ER7 = []; % {psi_0,psi_1}GRAD
MNPD = zeros(3600,1062);
MPD5 = MNPD; 
for ii=3:72
    load(strcat(d(ii).name,'/VertStruct.mat'),'NORM_PD','PD_COEF2',...
        'PD_COEF3','PD_COEF4')
     PD_COEF2 = PD_COEF2(:,:,2); % g/cm^3
%     PD_COEF4 = PD_COEF4(:,:,2); % g/cm^3
    load(strcat(d(ii).name,'/VertStruct2.mat'),'PD_COEF5')
    PD_COEF5 = PD_COEF5(:,:,2); % g/cm^3
    load(strcat(d(ii).name,'/Medley.mat'),'Err_P')
    
    MNPD = MNPD + NORM_PD;
    MPD5 = MPD5 + PD_COEF5;

    NORM_PD = mask(:,221:1000).*NORM_PD(:,221:1000);
    PD_COEF2 = mask(:,221:1000).*PD_COEF2(:,221:1000);
%     PD_COEF3 = mask(:,221:1000).*PD_COEF3(:,221:1000);
%     PD_COEF4 = mask(:,221:1000).*PD_COEF4(:,221:1000);
    PD_COEF5 = mask(:,221:1000).*PD_COEF5(:,221:1000);
    ERR_PD6 = mask(:,221:1000).*Err_P(:,221:1000,1);
    ERR_PD7 = mask(:,221:1000).*Err_P(:,221:1000,2);

    NPD = [NPD;NORM_PD(NORM_PD>TOL)];
    PD2 = [PD2;PD_COEF2(NORM_PD>TOL).^2];
%     PD3 = [PD3;PD_COEF3(NORM_PD>TOL).^2];
%     PD4 = [PD4;PD_COEF4(NORM_PD>TOL).^2];
    PD5 = [PD5;PD_COEF5(NORM_PD>TOL).^2];
    ER6 = [ER6;ERR_PD6(NORM_PD>TOL)];
    ER7 = [ER7;ERR_PD7(NORM_PD>TOL)];

    fprintf('ii = %2d \n',ii-2)
end
MNPD = MNPD/70;
MPD5 = MPD5/70;
clear NORM_PD PD_COEF* Err_P
