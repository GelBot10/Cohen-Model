%% H Model using Cohen2014
% This relize more on Fourier transforms than Fresnel propogation.
% Removed for loop for speed and possible GPU application.
% Added sparse conversion
% Does not save every step.
% Dynamic L0x, only 1 pz at a time.
% For running through all pz values.  Constants are initialized by the
% Initialization file depending on system you want to model.
% Pzs are then determined and the program is run as a for loop for each Pz
% value.  A PSF is calculated and saved for each Pz value then cleared
% before the next iteration.
%% Initialization
%% Physical properties
NA=0.14; %0.5; % Numerical Aperature
n=1; %1.33; % index of refraction in object space
FnOb=(2*tan(asin(NA/n)))^(-1); % F/# of the objective
lambda=780*10^(-9); %780*10^(-9); % wavelength of probe light in meters
k=2*pi*n/lambda; % wavenumber
alpha=asin(NA/n); % half angle of the NA
M=.66; %20; % magnification
fobj=40*10^(-3); % focal length
Dobj=fobj/FnOb; % Diameter of the objective
%tube lens
ftl=M*fobj; % from paper
Dtl=ftl/(M*FnOb);

eps0=8.854*10^(-12); %permittivity of free space
sol=299792458; %speed of light
%% sensor/mcirolens properties
% Pixel dimensions of the sensor
as=5250; %resampled; non resampled 5369
bs=7574; %resampled; non resampled 7728
% Pixel pitch
dpix=4.57*10^(-6);
% ML properties
NumPix=14; % number of pixels behind a micro lens in one dimension.
dml=dpix*NumPix; %62.5*10^(-6); %L2x/Nmlx; % diameter of the microlens
FnMl=2; % F/# of the microlenses
fml=dml*FnMl; % focal length determined by F/# of the micro lens
%% P space
L1px=dml/M;
Npx=11;
dpx=L1px/Npx; % voxel side length (m)

px=-L1px/2+dpx/2:dpx:L1px/2; % this is centered
py=px;
%% Saving things
% Excel file name
filename='I:\MATLAB\Dissertation_Code\Cohen Model PSFs\hHtNotesv9.xlsx';
[xstat,sheets,xlF]=xlsfinfo(filename);
[xln,xlt,xlr]=xlsread(filename,'A1:A100');
% Folder

mkdir(strcat('I:\MATLAB\Dissertation_Code\Cohen Model PSFs\',...
    datestr(datetime('now'),'yyyymmdd'),'-',num2str(NA),'-',num2str(Npx)))
foldername=strcat('I:\MATLAB\Dissertation_Code\Cohen Model PSFs\',...
    datestr(datetime('now'),'yyyymmdd'),'-',num2str(NA),'-',num2str(Npx),'\');
% strcat('hHtv9_',datestr(datetime('now'),'yyyymmdd_HHMMss'))
%%
Pz=PzValues(lambda, NA, NumPix, n, 1);
nz=size(Pz);
%%
% initialize Isf
Isf=sparse([]);
for ll=1:nz(2)
pz=Pz(ll);
[PX,PY,PZ]=meshgrid(px,py,pz);
[a,b,c]=size(PX);
%% x space sizes
%U0 plane
Nx=2^8-1; %2^11-1; %2^12-1; %floor(Npx*Nmlx); Npx*Nmlx; %2^11-1; 2^12+1; 2^13-3; % %2047;  
L0x=sqrt(lambda*abs(pz)*Nx); %Ls/M; %1.2*10^(-3);

dxx=L0x/Nx; %dml/Npx/M; %L0x/Nx;
xx=-L0x/2+dxx/2:dxx:L0x/2;
yy=xx;
[X,Y]=meshgrid(xx,yy);
[ax,bx]=size(X);

%U1 plane
L1x=lambda*fobj/dxx;
dxx1=lambda*fobj/L0x;
xx1=-L1x/2+dxx1/2:dxx1:L1x/2;
[X1,Y1]=meshgrid(xx1,xx1);

% U2 plane
L2x=lambda*ftl/dxx1;
dxx2=lambda*ftl/L1x;
xx2=-L2x/2+dxx2/2:dxx2:L2x/2;
[X2,Y2]=meshgrid(xx2,xx2);

%U3 plane, frequency domain only
fx=-1/(2*dxx2)+1/(2*L2x):1/L2x:1/(2*dxx2);
[FX,FY]=meshgrid(fx,fx);

%% Number of microlenses
Nmlx=floor((L0x*M)/dml); %floor(Ls/dml); Ls/dml; % number of mls allong a give direction
if mod(Nmlx,2)==0; Nmlx=Nmlx-1; end
Nmly=Nmlx; %floor(Ls/dml);%fax(2);
%% r
R=zeros(ax,bx,a,b,c);
for ii=1:a
    for jj=1:b
        for kk=1:c
            R(:,:,ii,jj,kk)=sqrt((X-PX(ii,jj,kk)).^2+(Y-PY(ii,jj,kk)).^2+PZ(ii,jj,kk)^2);
        end
    end
end

%% U_0
A=0.000001; % amplitude of point souce Units are such that the units of U0
% are V/m and watts per pixel detected is reasonable
U0=zeros(ax,bx,a,b,c);
for ii=1:c
    U0(:,:,:,:,ii)=-1i*A*k*PZ(1,1,ii)./(2*pi*(R(:,:,:,:,ii)).^2).*exp(1i*k*R(:,:,:,:,ii));%.*g0(:,:,ii,:,:);
end
clear R
%% Needed Definitions
% n=1; %now in air
% k=2*pi*n/lambda;


% Pre alocate U1m
% U1m=zeros(size(U0));
% Objective Pupil function
Pobj=rect(sqrt(X.^2+Y.^2)/(Dobj));
Pobj=repmat(Pobj,1,1,a,b,c);
% The Tube Lens and Pupil function
% Angular limit
% Ra1=sqrt(X1.^2+Y1.^2)/fobj;
% Pobja=zeros(size(Ra1));
% Pobja(Ra1<=alpha)=1;
 %fobj/(M/(2*NA));
Ptl=round(rect(sqrt(X1.^2+Y1.^2)/(Dtl)),10);
Ptl=repmat(Ptl,1,1,a,b,c);
% Pre alocate U2m
% U2m=zeros(size(U0));

% MLA phase
phi=rect(X2/dml).*rect(Y2/dml).*exp((-1i*k/(2*fml))*(X2.^2+Y2.^2)); %phase at the microlens
Npix=ceil(dml/dxx2);
[Gx,Gy]=meshgrid(1:ax,1:bx);
if mod(ax,2)==0
    cGx=rem(Gx-(ax)/2,Npix)==0;
    cGy=rem(Gy-(ax)/2,Npix)==0;
else
    cGx=rem(Gx-(ax+1)/2,Npix)==0;
    cGy=rem(Gy-(ax+1)/2,Npix)==0;
end
combX=cGx.*cGy; %ucomb(X2/dml).*ucomb(Y2/dml);
Phi=conv2(phi,combX,'same');
Phi=repmat(Phi,1,1,a,b,c);
% Pre alocate h
% h=zeros(a,b,c,ax,bx);
%% Prop

% for ii=1:a
%     for jj=1:b
%         for kk=1:c
            %Propogation through Objective
%             U0p=squeeze(U0(ii,jj,kk,:,:));%.*CosTh;%.*Pobj;
            
            ph=1/(1i*lambda*fobj)*exp(1i*k*fobj);

            U0=ph*fftshift(fftshift(fft2(ifftshift(ifftshift(U0,1),2)),1),2)*dxx^2;
            U0=U0.*Pobj;
            
            %Propogation through Tube lense
            U0=U0.*Ptl;
            ph=1/(1i*lambda*ftl)*exp(1i*k*ftl);
            
            U0=ph.*fftshift(fftshift(fft2(ifftshift(ifftshift(U0,1),2)),1),2)*dxx1^2;
            
            %Propogation through MLA to sensor
            U0=Phi.*U0;
            
            H=exp(-1i*pi*lambda*fml*(FX.^2+FY.^2));    %trans fucn
            H=repmat(H,1,1,a,b,c);
            H=ifftshift(ifftshift(H,1),2);  %shift trans fucn
            U0=fft2(U0);%fft2(fftshift(u1));  %shift, fft src field
            U0=H.*U0;   %multiply
            U0=ifft2(U0);%ifftshift(ifft2(U2));    %inv fft, center obs field
            
%         end
%     end
% end

% Ht=abs(U3).^2; %PSF
U0=eps0*sol*((abs(U0).^2)/2)*dxx2^2; %intensity in Watts per pixel
clear Phi H Pobj Ptl
% Sample down to actual pixel size
fun = @(block_struct) mean( block_struct.data(:) );
bsn=round(Npix/14); %Npixs);
ab=ceil(ax/bsn);
bb=ceil(bx/bsn);
if mod(ab,2)==0
    abm=ab/2;
    abl=as/2-abm+1;
    abh=as/2+abm;
else
    abm=(ab+1)/2;
    abl=as/2-abm+1;
    abh=as/2+abm-1;
end
if mod(bb,2)==0
    bbm=bb/2;
    bbl=bs/2-bbm+1;
    bbh=bs/2+bbm;
else
    bbm=(bb+1)/2;
    bbl=bs/2-bbm+1;
    bbh=bs/2+bbm-1;
end

Iz=zeros(as,bs,a,b,c);
for ii=1:a
    for jj=1:b
        for kk=1:c
            B=blockproc(U0(:,:,ii,jj,kk),[bsn bsn],fun);
            Iz(abl:abh,bbl:bbh,ii,jj,kk)=B;
        end
    end
end
% Iz=reshape(Iz,as*bs,a*b*c);
% Iz=sparse(Iz);
Imax=max(Iz(:)); %max intesity measured
Imin=Imax/255; % minimum intensity measurable
Isp=Iz; %reshape(Iz,as*bs,a*b*c);
Isp(Isp<Imin)=0; %set values under minimum to zero
Is=Isp; %sparse(Isp);
Isf=[Isf Is];
%%
label={strcat('U0-Iz-Is-v9_',num2str(ll))};
save(strcat(foldername,label{1}),'U0','Iz','Is','-v7.3')

end
%% Exporting
save(strcat(foldername,'Isf'),'Isf')
ExpVals=[NA,n,FnOb,lambda,M,fobj,L1px,Npx,L0x,Nx,Nmlx,Nmly,dml,fml,Npix,L2x,dxx2,A];
%'Format','yyyyMMdd_HHmmss'
label=strcat('IniVals-v9_',datestr(datetime('now'),'yyyymmdd_HHMMss'));
for ii=1:100
    if isnan(xlr{ii})
        xlswrite(filename,ExpVals,'Sheet1',strcat('B',num2str(ii)))
        xlswrite(filename,label,'Sheet1',strcat('A',num2str(ii)))
        xlswrite(filename,Pz,'Sheet2',strcat('B',num2str(ii)))
        xlswrite(filename,label,'Sheet2',strcat('A',num2str(ii)))
        break
    end
end