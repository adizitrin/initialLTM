function initial_LTM_model(member_catalog_file)

%This function is creating a so-called light-trcaes-mass mass model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written originally by Adi Zitrin and Tom Broadhurst (Zitrin et al., 2009, MNRAS, 396, 1985)
%Updated significantly since. Last updated: May 2019.
%Suggested and implemented improvements by other contributors is highly
%appreciated especially Matthias Bartelmann, Gregor Seidel, Irene Sendra, Keiichi Umetsu
%and much advice from CLASH, Hubble Frontier Fields, and RELICS collaborators.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters to be defined by User:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
No_of_Bright_Free=2; %no. of brightest galaxies with M/L weight to be modified
Ind_Free_Gal(1:No_of_Bright_Free)=[1 2]; %indices in the member catalog of thos brightest galaxies

%The following defin model's the input parameters (user can edit):
q=1.3; %power-law exponent
s=150; % msoothing length in pixels
k_new=1.5; %overall scaling
k_gal=0.07; %galaxy contiburion compared to dark matter
gamma=0; %external shear strength
phi=0; %external shear angle (in radians)
rci=0; %BCG core - assumes the BCG is the first line of member galaxy file
BCGel=0; %BCG ellipticity - can also be read from file (assumes the BCG is the first line of member galaxy file)
BCGPA=0; %BCG position angle in degrees - can also be read from file (assumes the BCG is the first line of member galaxy file)
rci2=0; %2nd BCG core - assumes the BCG is the second line of member galaxy file
BCGel2=0; %2nd BCG ellipticity - can also be read from file (assumes this BCG is the second line of member galaxy file)
BCGPA2=0; %2nd BCG position angle in degrees - can also be read from file (assumes this BCG is the second line of member galaxy file)

BCGweight=5;%A factor describing first BCG M/L ratio compared to all other galaxies
BCG2weight=3;%A factor describing second BCG M/L ratio compared to all other galaxies
%(recall to add more weight if more galaxies are left free)
%define (don't touch):
freegal_weight=[BCGweight BCG2weight];

pix_scale=0.065; % arcsec/pixel. Pixel scale for CLASH images used here
i_x=1000;% x pixel coordinate of start of FOV in the original HST image the photometry was made on).
i_y=1000;% y pixel coordinate of start of FOV in the original HST image the photometry was made on).
leng=3000;% Desired modelinf FOV size in HST pixels; NOTE: too small a FOV may and will cause artefacts

SplineOrGaussian='Gaussian'; %spline interpolation or gaussain smoothing
Wresolution=4;% a factor describing by how much the model resolution will initially be (later interpolated. Recommended to use 4 or lower here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some definitions (do not touch):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_data=load(member_catalog_file); %loading the member galaxy clatalog. The format of the member galaxy catalog should be:
%x and y pixels coordinates in columns 2 and 3
%flux or luminosity to be used, in column 7

resto=Wresolution;
x_size=i_x+leng;
y_size=i_y+leng;
leng2=round(leng/resto);

alpha_xf_ram(1:3*leng2, 1:3*leng2)=0;
[Y,X]=meshgrid(-leng:resto:leng-1,-leng:resto:leng-1);
[Yg,Xg]=meshgrid(-1.5*leng:resto:1.5*leng-1,-1.5*leng:resto:1.5*leng-1);
X1=(1:(x_size-i_x-1)/(resto)+1);
Y1=(1:(y_size-i_y-1)/(resto)+1);
[YI2,XI2]=meshgrid(1:resto:(y_size-i_y),1:resto:(x_size-i_x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
%Start calculations. First let's calculate the two BCGs deflection field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=1;
epsgal=BCGel;
PAgal=BCGPA*pi/180;
bcg_x=(file_data(t,2)-i_x)+1;
bcg_y=(file_data(t,3)-i_y)+1;
delta_x1=(XI2-bcg_x)*cos(PAgal)+(YI2-bcg_y)*sin(PAgal);
delta_y1=-(XI2-bcg_x)*sin(PAgal)+(YI2-bcg_y)*cos(PAgal);
rtemp=sqrt(1/(1+epsgal)^2*delta_x1.^2+1/(1-epsgal)^2*delta_y1.^2);
if find(Ind_Free_Gal(:)==t)>0
    atempto=find(Ind_Free_Gal(:)==t);
    sigmaa= file_data(t,7)*freegal_weight(atempto)*sqrt(rtemp.^2+rci.^2).^(-q);
else
    sigmaa= file_data(t,7)*sqrt(rtemp.^2+rci.^2).^(-q);
end

x_factor=X./(X.^2+Y.^2);
y_factor=Y./(X.^2+Y.^2);

x_factor_shift=fftshift(x_factor);
y_factor_shift=fftshift(y_factor);

x_factor_shift(1,1)=0;
y_factor_shift(1,1)=0;
xfourier=fft2(x_factor_shift);
yfourier=fft2(y_factor_shift);
smooth_comp=sigmaa;
new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;%mean2(smooth_comp_first);
new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;

fourier_smooth_new=fft2(new_smooth);

alpha_x_DMf=ifft2(xfourier.*fourier_smooth_new);
alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);

alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_x1=alpha_x_DMfs*resto^2/pi;
alpha_y1=alpha_y_DMfs*resto^2/pi;

%Second BCG:
t=2;
epsgal=BCGel2;
PAgal=BCGPA2*pi/180;
bcg_x2=(file_data(t,2)-i_x)+1;%file_data(t,2)-i_x+1;
bcg_y2=(file_data(t,3)-i_y)+1;%file_data(t,3)-i_y+1;
delta_x1=(XI2-bcg_x2)*cos(PAgal)+(YI2-bcg_y2)*sin(PAgal);%*cos(phi);
delta_y1=-(XI2-bcg_x2)*sin(PAgal)+(YI2-bcg_y2)*cos(PAgal);%*sin(phi);
rtemp=sqrt(1/(1+epsgal)^2*delta_x1.^2+1/(1-epsgal)^2*delta_y1.^2);
if find(Ind_Free_Gal(:)==t)>0
    atempto=find(Ind_Free_Gal(:)==t);
    sigmaa= file_data(t,7)*freegal_weight(atempto)*sqrt(rtemp.^2+rci2.^2).^(-q);
else
    sigmaa= file_data(t,7)*sqrt(rtemp.^2+rci2.^2).^(-q);
end

smooth_comp=sigmaa;
new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;
new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;

fourier_smooth_new=fft2(new_smooth);

alpha_x_DMf=ifft2(xfourier.*fourier_smooth_new);
alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);

alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
alpha_x2=alpha_x_DMfs*resto^2/pi;
alpha_y2=alpha_y_DMfs*resto^2/pi;

clear alpha_xf

alpha_xf=(alpha_xf_ram);

%Calculation of deflection field from rest of power-law galaxies
for t=3:length(file_data(:,1))
    if file_data(t,2)>i_x+1-leng2 && file_data(t,2)<i_x+leng-1+leng2 && file_data(t,3)>i_y+1-leng2 && file_data(t,3)<i_y+leng-1+leng2
        if find(Ind_Free_Gal(:)==t)>0
            atempto=find(Ind_Free_Gal(:)==t);
            alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7)*freegal_weight(atempto);
        else
            alpha_xf(round((file_data(t,2)-i_x+leng)/resto),round((file_data(t,3)-i_y+leng)/resto))=file_data(t,7);
            
        end
    end
end
alpha_yf=alpha_xf;

xf_factor=Xg./(Xg.^2+Yg.^2).^(0.5*q);
yf_factor=Yg./(Xg.^2+Yg.^2).^(0.5*q);

xf_factor_shift=fftshift(xf_factor);
yf_factor_shift=fftshift(yf_factor);

xf_factor_shift(1,1)=0;
yf_factor_shift(1,1)=0;
xffourier=fft2(xf_factor_shift);
yffourier=fft2(yf_factor_shift);
fourier_alpha_xf=fft2(alpha_xf);
fourier_alpha_yf=fft2(alpha_yf);

alpha_x_ff=ifft2(xffourier.*fourier_alpha_xf);
alpha_y_ff=ifft2(yffourier.*fourier_alpha_yf);

alpha_x_ffs=alpha_x_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2)*pi;
alpha_y_ffs=alpha_y_ff(leng2+1:leng2+leng2,leng2+1:leng2+leng2)*pi;
clear alpha_x_ff alpha_y_ff
alpha_x=real((alpha_x_ffs));
alpha_y=real((alpha_y_ffs));
alpha_x=alpha_x+alpha_x1+alpha_x2;
alpha_y=alpha_y+alpha_y1+alpha_y2;

%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the mass density distribution:
[da_x_dy,da_x_dx]=gradient(alpha_x/resto);
[da_y_dy,da_y_dx]=gradient(alpha_y/resto);
kappagals=(da_x_dx+da_y_dy)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In case wish to see mass distribution from galaxies uncomment:
%image(fliplr(poisson*2)')
%waitforbuttonpress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now smooth DM - with Spline interpolation or better off Gaussian kernel:
if strcmp(SplineOrGaussian,'Spline')==1
    W=kappagals;
    W(:,:)=1;
    ZI=polyfitweighted2(Y1,X1,kappagals,s,W);
    smooth_comp=polyval2(ZI,Y1,X1);
    [m,n]=size(kappagals);
    new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;%mean2(smooth_comp_first);
    new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;
    fourier_smooth_new=fft2(new_smooth);
    alpha_x_DMf=ifft2(xfourier.*fourier_smooth_new);
    alpha_y_DMf=ifft2(yfourier.*fourier_smooth_new);
    alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
    alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
    clear alpha_x_DMf alpha_y_DMf
    alpha_x_DM=alpha_x_DMfs*resto^2/pi;
    alpha_y_DM=alpha_y_DMfs*resto^2/pi;
    clear  alpha_x_DMfs alpha_y_DMfs
    
    
elseif strcmp(SplineOrGaussian,'Gaussian')==1
    smooth_comp=kappagals;
    
    new_smooth(1:2*length(smooth_comp(:,1)), 1:2*length(smooth_comp(1,:)))=0;
    new_smooth(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)))=smooth_comp;
    
    fourier_smooth_new=fft2(new_smooth);
    g_factor=exp(-(X.^2+Y.^2)./(2.0*s*s))./(2.0*pi*s*s);
    g_factor_shift=fftshift(g_factor);
    gfourier=fft2(g_factor_shift);
    alpha_x_DMf=ifft2(xfourier.*gfourier.*fourier_smooth_new);
    alpha_y_DMf=ifft2(yfourier.*gfourier.*fourier_smooth_new);
    
    alpha_x_DMfs=alpha_x_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
    alpha_y_DMfs=alpha_y_DMf(1:length(smooth_comp(:,1)),1:length(smooth_comp(1,:)));
    clear alpha_x_DMf alpha_y_DMf
    
    alpha_x_DM=alpha_x_DMfs*resto^2/pi;
    alpha_y_DM=alpha_y_DMfs*resto^2/pi;
    
    clear  alpha_x_DMfs alpha_y_DMfs
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set both deflection fields -- galaxies and DM -- to the same arbitrary
%scale of, on average, 200 pixels
rat=200/mean(mean(abs(alpha_x_DM)));
rat_gal=200/mean(mean(abs(alpha_x)));

%Calculate External shear:
gammacos2phi=gamma*cos(2*phi);
gammasin2phi=gamma*sin(2*phi);
alpha_x_add=gammacos2phi*XI2+gammasin2phi*YI2;
alpha_y_add=gammasin2phi*XI2-gammacos2phi*YI2;

%Interpolate to original resolution:
clear alpha_x_ALL1 alpha_y_ALL1
alpha_x_ALL1=(k_gal*rat_gal*alpha_x(1:leng2,1:leng2)+(1-k_gal)*rat*alpha_x_DM(1:leng2,1:leng2))*k_new+alpha_x_add(1:leng2,1:leng2);
alpha_y_ALL1=(k_gal*rat_gal*alpha_y(1:leng2,1:leng2)+(1-k_gal)*rat*alpha_y_DM(1:leng2,1:leng2))*k_new+alpha_y_add(1:leng2,1:leng2);
alpha_x_ALL2=interp2(alpha_x_ALL1,1:(1/Wresolution):leng2+(1-1/Wresolution),(1:(1/Wresolution):leng2+(1-1/Wresolution))');
alpha_y_ALL2=interp2(alpha_y_ALL1,1:(1/Wresolution):leng2+(1-1/Wresolution),(1:(1/Wresolution):leng2+(1-1/Wresolution))');

alpha_x_ALL=(alpha_x_ALL2);
alpha_y_ALL=(alpha_y_ALL2);
clear alpha_x_ALL1 alpha_y_ALL1 alpha_x_ALL2 alpha_y_ALL2

%Get mass distribtuion and magnification:
[da_x_dy,da_x_dx]=gradient(alpha_x_ALL);
[da_y_dy,da_y_dx]=gradient(alpha_y_ALL);
poisson_ALL=da_x_dx+da_y_dy;
magnification_ALL=abs(1./(1-poisson_ALL+da_x_dx.*da_y_dy-da_x_dy.*da_y_dx));
close all

%Show the kappa maps and magnification map:
figure;
subplot(1,2,1)
image(fliplr((poisson_ALL/2)*40)'); title('kappa map'); %the *40 is arbitrary factor for showing image purposes
subplot(1,2,2)
image(fliplr((magnification_ALL))'); title('magnification map')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('initialLTMmodel.mat','alpha_x_ALL','alpha_y_ALL','i_x','i_y','leng','x_size','y_size');