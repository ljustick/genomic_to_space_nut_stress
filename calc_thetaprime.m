function x=calc_thetaprime(chl,bbp,par,mld,yd)

%%% matlab fxn to calculate theta' given inputs of Chl, bbp(443), PAR, MLD, and day-of-year (1-365).
%%% can work on individual points or on matrices (e.g., satellite fields).  For our
%%% efforts, we used input fields of NASA Level 3 mapped products (2160x4320)
%%% equivalent to ~9.25km resolution at the equator.


%%% -- define 9km lat/lon grid  and calculate daylength based on lat and day-of-year
%%% -- daylength sub-function appended below

lat9=90-1/24:-1/12:-90+1/24; 
lon9=-180+1/24:1/12:180-1/24;         %%% -- ~9km res

[lon9g,lat9g]=meshgrid(lon9,lat9);
dl=daylength(yd,lat9g);

%%% -- calcluate Cphyto from bbp(443) following Behrenfeld et al. (2005) and resultant C:Chl ratio
Cphy=(bbp-0.00035).*13000; 
Cphy(Cphy<0)=NaN;
cc=Cphy./chl;

%%% -- calculate average attenuation of PAR in euphotic zone (Kpar)
%%% -- from Morel et al. (2007)
c=log10(chl);
logzeu=1.524 - 0.436.*c - 0.0145.*(c.^2) + 0.0186.*(c.^3); 
zeu=10.^logzeu;
kpar=-log(0.01)./zeu;

%%% -- calculate median mixed layer growth irradiace (Ig) and photoacclimation-based C:Chl following
%%% -- Behrenfeld et al. (2015)
parDL=par./dl;
Ig=0.98.*parDL.*exp(-kpar.*(mld./2));
chlCSM=(1+exp(-0.15.*parDL))./(1+exp(-3.*Ig));
CChlmod=(19.*exp(0.038.*(parDL.^0.45)./kpar)).*chlCSM;

%%% -- thetaprime = C:Chl normalized to photoacclimation-based C:Chl
x=cc./CChlmod;

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