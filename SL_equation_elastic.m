% Code to solve the elastic sea level equation following 
% Kendall et al., 2005 and Austermann et al., 2015

% J. Austermann 2015

% add paths when run for the first time.
% addpath SLFunctions
% addpath '/Users/jackyaustermann/Documents/MATLAB/m_map'

%% Parameters & Input 
% Specify maximum degree to which spherical transformations should be done
maxdeg = 512;

% parameters
rho_ice = 916;
rho_water = 1000;
rho_sed = 2300;
g = 9.81;


% The following steps help speed up the calculations
% Set up Gauss Legendre grid onto which to interpolate all grids
N = maxdeg; 
[x,w] = GaussQuad(N);
x_GL = acos(x)*180/pi - 90;
lon_GL = linspace(0,360,2*N+1);
lon_GL = lon_GL(1:end-1);

colat = 90 - x_GL;
lon = lon_GL;

[lon_out,lat_out] = meshgrid(lon_GL,x_GL);

% Precompute legendre polynomials
P_lm = cell(N+1,1);
for l=0:N
    P_lm{l+1} = legendre(l,x,'norm');
end


% --------------------------------
% ICE
% --------------------------------

load WAIS 

ice_0_nointerp = ice_Ant;
ice_j_nointerp = ice_EAIS;

% interpolate ice masks on common grid
ice_0 = interp2(lon_WAIS,lat_WAIS,ice_0_nointerp,lon_out, lat_out);
ice_j = interp2(lon_WAIS,lat_WAIS,ice_j_nointerp,lon_out, lat_out);

del_ice = ice_j - ice_0; 


% --------------------------------
% DYNAMIC TOPOGRAPHY
% --------------------------------

del_DT = zeros(size(del_ice));


% --------------------------------
% SEDIMENT
% --------------------------------

del_sed = zeros(size(del_ice));


% --------------------------------
% TOPOGRAPHY
% --------------------------------

% load preloaded etopo (including ice) as topo_orig, lon_topo, lat_topo
load topo_SL

% interpolate topography grid onto Gauss Legendre Grid
topo0 = interp2(lon_topo,lat_topo,topo_orig,lon_out, lat_out);



%% Set up love number input

% prepare love numbers in suitable format and calculate T_lm and E_lm 
% to calculate the fluid case, switch h_el to h_fl, k_el to k_fl and same
% for tidal love numbers
load SavedLN/LN_l90_VM2
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
h_lm_tide = love_lm(h_el_tide,maxdeg);
k_lm_tide = love_lm(k_el_tide,maxdeg);

E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);

E_lm_T = 1 + k_lm_tide - h_lm_tide;

% can switch this in if you want to exclude rotational effects
% E_lm_T = zeros(size(E_lm_T));

%% Solve sea level equation (after Kendall 2005, Dalca 2013 & Austermann et al. 2015)

k_max = 10;   % maximum number of iterations
epsilon = 10^-4; % convergence criterion

% 0 = before
% j = after

% set up present-day topo and ocean function 
topo_0 = topo0; % already includes ice and dynamic topography
oc_0 = sign_01(topo_0);

% set up topography and ocean function after the ice change
topo_j = topo_0 + del_ice; % del_ice is negative -> subtract ice that is melted
oc_j = sign_01(topo_j);

% calculate change in sediments and decompose into spherical harmonics
Sed_lm = sphere_har(del_sed,maxdeg,N,P_lm); 

% expand ocean function into spherical harmonics
oc0_lm = sphere_har(oc_0,maxdeg,N,P_lm); 


% initial values for convergence
conv = 'not converged yet';

        
for k = 1:k_max % loop for sea level and topography iteration

    switch conv

        case 'converged!'

        case 'not converged yet'
            
        % expand ocean function into spherical harmonics
        ocj_lm = sphere_har(oc_j,maxdeg,N,P_lm);  
        
        % CHECK ICE MODEL 
        % check ice model for floating ice
%         check1 = sign_01(-topo_j + ice_j);
%         check2 = sign_01(+topo_j - ice_j) .* ...
%          (sign_01(-ice_j*rho_ice - (topo_j - ice_j)*rho_water));
%         
%         ice_j_corr = check1.*ice_j + check2.*ice_j;
%         del_ice_corrected = ice_j_corr - ice_0; 
        del_ice_corrected = ice_j - ice_0; 
        
        deli_lm = sphere_har(del_ice_corrected,maxdeg,N,P_lm);  
        
        
        % calculate topography correction
        TO = topo_0.*(oc_j-oc_0);
        % expand TO function into spherical harmonics
        TO_lm = sphere_har(TO,maxdeg,N,P_lm); 
        
        
        % set up initial guess for sea level change
        if k == 1
            % initial guess of sea level change is just to distribute the
            % ice over the oceans
            delS_lm = ocj_lm/ocj_lm(1)*(-rho_ice/rho_water*deli_lm(1) + ...
                TO_lm(1));
            % convert into spherical harmonics
            % delS_init = inv_sphere_har(delS_lm,maxdeg,N,P_lm);
            
        end
        
        % calculate loading term
        L_lm = rho_ice*deli_lm + rho_water*delS_lm + rho_sed*Sed_lm;

        % calculate contribution from rotation
        La_lm = calc_rot(L_lm,k_el,k_el_tide);

        % calculate sea level perturbation
        % add ice and sea level and multiply with love numbers
        % DT doesn't load!
        delSLcurl_lm_fl = E_lm .* T_lm .* (rho_ice*deli_lm + rho_water*delS_lm + rho_sed*Sed_lm) + ...
            1/g*E_lm_T.*La_lm;

        % convert to spherical harmonics and subtract terms that are part
        % of the topography to get the 'pure' sea level change
        delSLcurl_fl = sph2spa(delSLcurl_lm_fl,maxdeg,lon,colat,P_lm_sph2spa);
        delSLcurl = delSLcurl_fl - del_ice_corrected - del_DT - del_sed;


        % compute and decompose RO
        RO = delSLcurl.*oc_j;
        RO_lm = spa2sph(RO,maxdeg,lon,colat,P_lm_spa2sph);

        % calculate eustatic sea level perturbation (delta Phi / g)
        delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*deli_lm(1) ...
            - RO_lm(1) + TO_lm(1));

        % calculate overall perturbation of sea level over oceans
        % (spatially varying field and constant offset)
        delSL = delSLcurl + delPhi_g;


        % update topography and ocean function
        topo_j = - delSL + topo_0;
        oc_j = sign_01(topo_j);


        % calculate change in ocean height and decompose
        delS_new = delSL.*oc_j -  topo_0.*(oc_j-oc_0);
        delS_lm_new = spa2sph(delS_new,maxdeg,lon,colat,P_lm_spa2sph);


        % calculate convergence criterion chi
        chi = abs( (sum(abs(delS_lm_new)) - sum(abs(delS_lm))) / ...
            sum(abs(delS_lm)) );

        % check convergence against the value epsilon
        % If converged, set the variable conv to 'converged!' so that the
        % calculation exits the loop. If not converged iterate again.
        if chi < epsilon;
            conv = 'converged!';
            disp(['Converged after iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
        else
            conv = 'not converged yet';
            disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
        end

        % update sea sea surface height
        delS_lm = delS_lm_new;
    end

end

% calculate the scaling to normalize the fingerprint (it's normalized to be
% one on average, when averaged over the final ocean basin). 
% calculate change in sea level over final ocean basin
del_scaling =(delSL + del_ice_corrected).*oc_j;
% get the average of that when spreading the water over the whole globe
sca = spa2sph(del_scaling,1,lon,colat);
% get the average of that when spreading the water only over the oceans.
scaling_fact = sca(1)/ocj_lm(1);  


%% Plot results

% We only want the sea level change cause by melted ice, so subtract
% del_ice
SL_change_rot = delSL + del_ice_corrected;
% normalize it by the scaling factor
plotSL = SL_change_rot/scaling_fact;

% construct identical colormap to Mitrovica 2009 paper
MyColorMap = [238   44  37
    211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244
    173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235
    163 201 235
    111 147 201
    96  103 175
    74  102 176
    68  87  165
    58  84  163
    53  69  154
    44  47  137
    38  35  103
    19  15  54
    0   0   0
    0   0   0];

% plot
figure
m_proj('hammer-aitoff','clongitude',0);
m_pcolor([lon_out(:,end/2+1:end)-360 lon_out(:,1:end/2)],lat_out,...
    [plotSL(:,end/2+1:end) plotSL(:,1:end/2)])
m_coast('color',[0 0 0]);
m_grid('box','fancy','xticklabels',[],'yticklabels',[]);
shading flat
colorbar
colormap(MyColorMap/255)
caxis( [-0.05, 1.6] )
title('Sea level fingerprint of West Antarctic Ice Sheet collapse')

