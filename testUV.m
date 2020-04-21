load Depth_Land_Masks.mat mask % the mask array has 1 in locations with > 2km depth and >200 points from shore
cd 1Deg/
d = dir;
cd ../
TOL = 4; % (cm/s)^2
EKE = [];
ERR1 = [];
ERR2 = [];
ERR3 = [];
ERR4 = [];
MEKE = zeros(3600,1062);
MER1 = MEKE;
for ii=3:72
    load(strcat(d(ii).name,'/VertStruct.mat'),'NORM_UV')
    load(strcat(d(ii).name,'/Kappa_EKE.mat'),'Err_K')
    tmp = Err_K;
    load(strcat(d(ii).name,'/Medley.mat'),'Err_K')
    Err_K = cat(3,tmp,Err_K);clear tmp

    MEKE = MEKE + NORM_UV;
    MER1 = MER1 + Err_K(:,:,1);

    NORM_UV = mask(:,221:1000).*NORM_UV(:,221:1000);
    Err_K = bsxfun(@times,mask(:,221:1000),Err_K(:,221:1000,:));

    Err_1 = Err_K(:,:,1);
    Err_2 = Err_K(:,:,2);
    Err_3 = Err_K(:,:,3);
    Err_4 = Err_K(:,:,4);clear Err_K

    EKE = [EKE;NORM_UV(NORM_UV>TOL)];
    ERR1 = [ERR1;Err_1(NORM_UV>TOL)];
    ERR2 = [ERR2;Err_2(NORM_UV>TOL)];
    ERR3 = [ERR3;Err_3(NORM_UV>TOL)];
    ERR4 = [ERR4;Err_4(NORM_UV>TOL)];

    fprintf('ii = %2d \n',ii-2)
end
MEKE = MEKE/70;
MER1 = MER1/70;
clear NORM_UV Err_*
