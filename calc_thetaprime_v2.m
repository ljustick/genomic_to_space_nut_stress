function x=calc_thetaprime(chl,bbp,k490,par,mld,yd);

%%% matlab fxn to calculate theta' given inputs of Chl, bbp(532), k490, PAR, MLD, and day-of-year (1-365).
%%% This was designed to work on NASA Level 3 mapped products (2160x4320) and comparable MLD data field
%%%
%%%

%%% -- define 9km lat/lon grid and calculate daylength based on lat and day-of-year
%%% -- daylength sub-function appended below

lat9=90-1/24:-1/12:-90+1/24;, lon9=-180+1/24:1/12:180-1/24;         %%% -- ~9.25km res
[lon9g,lat9g]=meshgrid(lon9,lat9);
dl=daylength(yd,lat9g);

%%% -- calcluate Cphyto from bbp following Behrenfeld et al. (2023). Note that the input is bbp(532)
%%% -- and can be calculated using the bbp(443) and bbp spectral slope products and a typical power law
%%% -- this is also why the y-intercept (=0.00029) is different from Behrenfeld et al. (2005)
Cphy=(bbp-0.00029).*13000;, Cphy(Cphy<0)=NaN;
CChl_obs=Cphy./chl; %%  resultant C:Chl ratio

%%% -- calculate average attenuation of PAR in euphotic zone (Kpar)
%%% -- from Morel et al. (2007)
kpar=0.0665 + 0.874.*k490 - 0.00121./k490;

%%% -- calculate median mixed layer growth irradiace (Ig) and photoacclimation-based C:Chl 
%%% -- following Behrenfeld et al. (2015). A difference is that the coefficient in the deep-mixing
%%% -- solution (CChlDM) was 're-tuned' because of change from using GSM bbp in original publication 
%%% -- to GIOP bbp in current version. Exponent was changed to 40, then calculated  globally for the 
%%% -- whole MODIS time-series and then CChlDM was set to the average value (= 150).
parDL=par./dl;
Ig=0.98.*parDL.*exp(-kpar.*(mld./2));
CChlSM=(1+exp(-0.15.*parDL))./(1+exp(-3.*Ig));
CChlDM=150; %%% normally calculate as (40.*exp(0.038.*(parDL.^0.45)./kpar)), but 150 is the global average
CChl_photo=CChlDM.*chlCSM;

%%% -- thetaprime = (C:Chl)_photo : (C:Chl)obs
x=CChl_photo./CChl_obs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dl=daylength(yd,lat)
%
% calculates daylength.  required inputs are day of year (yd)
% and latitude (lat)
%
%
dec=23.5.*cos((2.*pi.*(yd - 172))./365);	%%%declination of sun
dl=acosd(-tand(lat).*tand(dec))./15.*2;		%%%sunrise hour angle / 15deg/hr x 2 (sr +ss)
dl=real(dl);                                    %%% real is in case 90 or -90 latitude is used
return