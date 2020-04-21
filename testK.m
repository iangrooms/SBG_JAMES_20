load Depth_Land_Masks.mat mask % the mask array has 1 in locations with > 2km depth and >200 points from shore
cd 1Deg/
d = dir;
cd ../
TOL = 4E-8; % (g/s cm^2)^2
NEF = [];
EK1 = [];
EK2 = []; % Aniso, k~1
EK3 = [];
EK4 = []; % Aniso, k~phi
EK5 = [];
EK6 = []; % Aniso, k~N^2
EK7 = []; % Ansio, k~EBT+
EK8 = []; % Aniso, k~BT+BC
for ii=3:72
    load(strcat(d(ii).name,'/VertStruct.mat'),'NORM_FR')
    load(strcat(d(ii).name,'/Kappa_EKE.mat'),'Err_F')
    load(strcat(d(ii).name,'/Medley.mat'),'Err_FK')

    NORM_FR = mask(:,221:1000).*NORM_FR(:,221:1000);
    Err1 = mask(:,221:1000).*Err_F(:,221:1000,1);
    Err2 = mask(:,221:1000).*Err_F(:,221:1000,2);
    Err3 = mask(:,221:1000).*Err_F(:,221:1000,3);
    Err4 = mask(:,221:1000).*Err_F(:,221:1000,4);
    Err5 = mask(:,221:1000).*Err_F(:,221:1000,5);
    Err6 = mask(:,221:1000).*Err_F(:,221:1000,6);
    Err7 = mask(:,221:1000).*Err_FK(:,221:1000,1);
    Err8 = mask(:,221:1000).*Err_FK(:,221:1000,2);
    
    NEF = [NEF;NORM_FR(NORM_FR>TOL)];
    EK1 = [EK1;Err1(NORM_FR>TOL)];
    EK2 = [EK2;Err2(NORM_FR>TOL)];
    EK3 = [EK3;Err3(NORM_FR>TOL)];
    EK4 = [EK4;Err4(NORM_FR>TOL)];
    EK5 = [EK5;Err5(NORM_FR>TOL)];
    EK6 = [EK6;Err6(NORM_FR>TOL)];
    EK7 = [EK7;Err7(NORM_FR>TOL)];
    EK8 = [EK8;Err8(NORM_FR>TOL)];

    fprintf('ii = %2d \n',ii-2)
end
clear NORM_FR K_X K_Y K_C Err*
