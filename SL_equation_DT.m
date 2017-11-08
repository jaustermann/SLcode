
%
%  Copyright (C) 2016 - 2017 by J. Austermann.
%  This file is part of SLcode.
%  SLcode is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%  SLcode is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  <http://www.gnu.org/licenses/>.

% Code to solve the sea level equation. Follows class notes of Jerry
% Mitrovica (lecture 17 of sea level notes) with some changes described by
% Dalca et al., 2013. The addition here is Dynamic topography
% Code solves the elastic/fluid case neglecting rotational effects and
% viscoelastic behaviour.
% J. Austermann 2013

% It is important to note that the RESULTS ARE IN 3 MA REFERENCE FRAME!
% 0 refers to present-day and j refers to a time in the past, e.g. 3 Ma

%% -----------------
% PARAMETERS & INPUT
% ------------------

function T_1 = SL_equation_DT(del_DT,del_S,del_I, ...
    del_G,T_0, maxdegree)

% Specify maximum degree to which spherical transformations should be done
maxdeg = maxdegree;

% parameters
rho_ice = 916;
rho_water = 1000;
rho_sed = 2300;
rho_air_mantle = 3300;
rho_water_mantle = 2300;
g = 9.81;

% Set up Gauss Legendre grid
N = maxdegree;
[x,w] = GaussQuad(N);
x_GL = acos(x)*180/pi - 90;
lon_GL = linspace(0,360,2*N+1);
lon_GL = lon_GL(1:end-1);

colat = 90 - x_GL;
lon = lon_GL;

[lon_out,lat_out] = meshgrid(lon_GL,x_GL);


%% -----------------------
% Set up love number input
% ------------------------

% load Love numbers
% prepare love numbers in suitable format and calculate T_lm and E_lm 
% load ReadSaveLN/LN_l70_ump5_lm5
load SavedLN/prem.l90C.umVM2.lmVM2.mat
h_lm = love_lm(h_fl, maxdeg);
k_lm = love_lm(k_fl, maxdeg);
h_lm_tide = love_lm(h_fl_tide,maxdeg);
k_lm_tide = love_lm(k_fl_tide,maxdeg);

E_lm = 1 + k_lm - h_lm;
T_lm = get_tlm(maxdeg);

E_lm_T = 1 + k_lm_tide - h_lm_tide;


%% ---------------------------------------------------------
% Solve sea level equation (after Kendall 2005 & Dalca 2013)
% ----------------------------------------------------------

% Calculation is in the 3 Ma reference frame. For SL change since 3 Ma need
% to use - delS etc. 

k_max = 10;         % maximum number of iterations
epsilon = 1e-4;    % convergence criterion


% set up present-day topo and ocean function (based on 3 Myr topo)
topo_0 = T_0;   % already includes ice, sediments and dynamic topography
oc_0 = sign_01(topo_0);

% set up topography and ocean function at later point in time
topo_j = T_0 + del_DT + del_S + del_I - del_G;
oc_j = sign_01(topo_j);


% Decompose change in sediments into spherical harmonics
Sed_lm = spa2sph(del_S,maxdeg,lon,colat);

% Expands grav_MC function into spherical harmonics
G_lm = spa2sph(del_G,maxdeg,lon,colat);

%ice_j_corr = check1.*ice_j + check2.*ice_j;
%delIce = ice_j_corr - ice_0;

delIce = del_I;

% Decompose change in ice
i_lm = spa2sph(delIce,maxdeg,lon,colat);

% Initial values for convergence
conv = 'not converged yet';
        
for k = 1:k_max % loop for sea level and topography iteration

    switch conv

        case 'converged!'

        case 'not converged yet'

        % Expand ocean function into spherical harmonics
        ocj_lm = spa2sph(oc_j,maxdeg,lon,colat);

        % Check ice model for floating ice
        %check1 = sign_01(-topo_j + ice_j);
        %check2 = sign_01(+topo_j - ice_j) .* ...
        % (sign_01(-ice_j*rho_ice - (topo_j - ice_j)*rho_water));
        
        %ice_j_corr = check1.*ice_j + check2.*ice_j;
        %delIce = ice_j_corr - ice_0;
        
        %delIce = del_I;
        
        % Decompose change in ice
        % i_lm = spa2sph(delIce,maxdeg,lon,colat_inv);
        
        % Calculate topography correction for this step
        topo_corr_j = topo_0.*(oc_j-oc_0);
        % Expand topo function into spherical harmonics
        TO_lm = spa2sph(topo_corr_j,maxdeg,lon,colat);


        if k == 1
            % Initial guess for change in sea surface height
            delSLcurl = topo_0 - topo_j; %= del_G - (del_H+del_I+del_DT)
            RO = delSLcurl.*oc_j;
            RO_lm = spa2sph(RO,maxdeg,lon,colat);
            
            delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*i_lm(1) ...
            - RO_lm(1) + TO_lm(1));
        
            delS_lm = RO_lm + delPhi_g * ocj_lm - TO_lm;
        end

        % Loading due to rotational feedback
        L_lm = rho_ice*i_lm + rho_water*delS_lm + rho_sed*Sed_lm;
        % Calculate effect of rotational change due to loading
        La_lm = calc_rot(L_lm,k_el,k_el_tide);
        
        % Calculate sea level perturbation
        % Add ice and sea level and multiply with love numbers
        % Don't include dynamic topography, because it doesn't load!
        delSLcurl_lm_fl = E_lm .* T_lm .* ...
            (rho_ice*i_lm + rho_water*delS_lm + rho_sed*Sed_lm) + ...
            1/g*E_lm_T.*La_lm + G_lm;
        delSLcurl_fl = sph2spa(delSLcurl_lm_fl,maxdeg,lon,colat);
        delSLcurl = delSLcurl_fl - del_I - del_DT - del_S;

        % Compute and decompose RO
        RO = delSLcurl.*oc_j;
        RO_lm = spa2sph(RO,maxdeg,lon,colat);

        % Compute delta Phi / g
        delPhi_g = 1/ocj_lm(1) * (- rho_ice/rho_water*i_lm(1) ...
            - RO_lm(1) + TO_lm(1));

        % Calculate overall sea level change (spatially varying and
        % constant component)
        delSL = delSLcurl + delPhi_g;

        % update topography and ocean function (sea level is the negative
        % of topography)
        topo_j =  topo_0 - delSL;
        oc_j = sign_01(topo_j);

        %compute and decompose ocean height
        delS_new = delSL.*oc_j -  topo_0.*(oc_j-oc_0);
        delS_lm_new = spa2sph(delS_new,maxdeg,lon,colat);

        % check convergence
        chi = abs( (sum(abs(delS_lm_new)) - sum(abs(delS_lm))) / ...
            sum(abs(delS_lm)) );

        if chi < epsilon;
            conv = 'converged!';
            disp(['Converged after iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
        else
            conv = 'not converged yet';
            disp(['Finished iteration ' num2str(k) '. Chi was ' num2str(chi) '.'])
        end

        % new ocean height for next iteration
        delS_lm = delS_lm_new;
    end

end

T_1 = T_0 - delSL;
end