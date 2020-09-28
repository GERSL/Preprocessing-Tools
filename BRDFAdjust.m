function ref_norm = BRDFAdjust( ref_ori,band_name,...
    solar_zenith_angle,view_zenith_angle,solar_azimuth_angle,view_azimuth_angle,...
    varargin)
% BRDFADJUST is to use the c-factor approach (Roy, D. P. et al., 2016) based on the RossThick-LiSparse-R BRDF model (Schaaf, Crystal B., et al. 2002).
%--------------------------------------------------------------------------
% Inputs:
% ref_ori: orginal surface reflectance
% band_name:  'Blue', 'Green', 'Red', 'NIR', 'SWIR1',  'SWIR2'
% solar_zenith_angle (unit: decimal degrees)
% view_zenith_angle (unit: decimal degrees)
% solar_azimuth_angle (unit: decimal degrees)
% view_azimuth_angle (unit: decimal degrees)
% centre_lat (unit: decimal degrees)
% solar_zenith_angle_norm (unit: decimal degrees)
%--------------------------------------------------------------------------
% Summary:
% This function is to adjust the BRDFs of Landsat and Sentinel images to
% niar observations.
%--------------------------------------------------------------------------
% History of revisions:
% use varargin by Shi Qiu 9/28/2020
% modify this function. by Rong Shang and Shi Qiu 6/27/2018
% create this function. by Shi Qiu and Zhanmang Liao 5/13/2017
% lamda=3;%Ñband number
%==========================================================================
    solar_zenith_angle = double(solar_zenith_angle);
    view_zenith_angle = double(view_zenith_angle);
    solar_azimuth_angle = double(solar_azimuth_angle);
    view_azimuth_angle = double(view_azimuth_angle);

    solar_zenith_angle = deg2rad(solar_zenith_angle);
    view_zenith_angle = deg2rad(view_zenith_angle);
    relative_azimuth_angle = deg2rad(view_azimuth_angle-solar_azimuth_angle);
    
    % if have centre's latitude
    p = inputParser;
    p.FunctionName = 'paras';
    % optional
    % default values.
    addParameter(p,'centre_lat', []);
    addParameter(p,'solar_zenith_angle_norm',[]);
    parse(p,varargin{:});
    centre_lat=p.Results.centre_lat;
    solar_zenith_angle_norm=p.Results.solar_zenith_angle_norm;
    if isempty(centre_lat)&&isempty(solar_zenith_angle_norm)
        error('No centre latitude or normalized solar zenith angle for BRDF correction\n');
    end
    
    if exist('centre_lat','var')&&~isempty(centre_lat)
        solar_zenith_angle_norm_centre = 31.0076 - 0.1272*centre_lat + 0.01187*centre_lat^2 + ...
            2.4*10^(-5)*centre_lat^3  - 9.48*10^(-7)*centre_lat^4 - ...
            1.95*10^(-9)*centre_lat^5 + 6.15*10^(-11)*centre_lat^6;
        solar_zenith_angle_norm_centre = deg2rad(solar_zenith_angle_norm_centre);
    end
    
    % if denfine solar angle
    if exist('solar_zenith_angle_norm','var')&&~isempty(solar_zenith_angle_norm)
        solar_zenith_angle_norm_custom = deg2rad(solar_zenith_angle_norm);
    end
    
    if exist('solar_zenith_angle_norm_custom','var')&&~isempty(solar_zenith_angle_norm_custom)
        % do not change the solar zenith angle if no any inputs because the
        % Landsat solar zenith angles are very similar.
        solar_zenith_angle_norm = solar_zenith_angle_norm_custom;
    else
        if exist('solar_zenith_angle_norm_centre','var')&&~isempty(solar_zenith_angle_norm_centre)
            solar_zenith_angle_norm = solar_zenith_angle_norm_centre;
        else
            % do not change solar angle if no self-defined one
            solar_zenith_angle_norm = solar_zenith_angle;
        end
    end
    
    view_zenith_angle_norm = deg2rad(0);
    relative_azimuth_angle_norm = deg2rad(180);

    % band name should be converted to band 1 2 3 4 5 and 7 for Landsat 7 data.
    switch band_name
        case 'Blue'
            lamda = 1; % band 1
        case 'Green'
            lamda = 2; % band 2
        case 'Red'
            lamda = 3; % band 3
        case 'NIR'
            lamda = 4; % band 4
        case 'SWIR1'
            lamda = 5; % band 5
        case 'SWIR2'
            lamda = 6; % band7
        otherwise
            ref_norm = [];
            warning('Please input right band name. (See BRDFAdjust)');
            return;
    end

    %%======== BRDF parameters  ========
    %  band   1       2       3       4       5       7
    f_iso= [0.0774, 0.1306, 0.1690, 0.3093, 0.3430, 0.2658];
    f_vol= [0.0372, 0.0580, 0.0574, 0.1535, 0.1154, 0.0639];
    f_geo= [0.0079, 0.0178, 0.0227, 0.0330, 0.0453, 0.0387];
    tmp = f_iso(lamda); clear f_iso;
    f_iso = tmp;
    tmp = f_geo(lamda); clear f_geo;
    f_geo = tmp;
    tmp = f_vol(lamda); clear f_vol;
    f_vol = tmp;
    %clear tmp;

    %======== calculate the kernel ========
    [k_vol_norm,k_geo_norm]=kernel(solar_zenith_angle_norm, view_zenith_angle_norm, relative_azimuth_angle_norm);
  %  clear solar_zenith_angle_norm view_zenith_angle_norm relative_azimuth_angle_norm;
    [k_vol_sensor,k_geo_sensor]=kernel(solar_zenith_angle, view_zenith_angle, relative_azimuth_angle);
 %   clear solar_zenith_angle view_zenith_angle relative_azimuth_angle;

    %======== calculate correcting parameter  ========
    P1=f_iso+f_geo*k_geo_norm+f_vol*k_vol_norm;
%    clear k_vol_norm;
    P2=f_iso+f_geo*k_geo_sensor+f_vol*k_vol_sensor;
%    clear k_vol_sensor;
    C_lamda=P1./P2;

    ref_norm = int16(C_lamda.*double(ref_ori)); %corrected reflectance

end

function [k_vol,k_geo] = kernel(theta_i,theta_v,azimuth_angle)
b=1;r=1;h=2;

%================ calculate k_vol  ================
cos_g = cos(theta_i).*cos(theta_v)+sin(theta_i).*sin(theta_v).*cos(azimuth_angle);
g = acos(max(-1,min(cos_g,1)));
clear cos_g;

k_vol = ((0.5*pi-g).*cos(g)+sin(g))./(cos(theta_i)+cos(theta_v))-pi/4;
clear g;


%================ calculate k_geo  ================
theta_i1 = atan(max(b/r*tan(theta_i),0));
theta_v1 = atan(max(b/r*tan(theta_v),0));
clear theta_i theta_v r;

g_1 = cos(theta_i1).*cos(theta_v1) + sin(theta_i1).*sin(theta_v1).*cos(azimuth_angle);
g_1 = acos( max(-1,min(1,g_1)) );

D = tan(theta_i1).^2+tan(theta_v1).^2-2*tan(theta_i1).*tan(theta_v1).*cos(azimuth_angle);
D = sqrt(max(D,0));

cos_t = h/b*( sqrt(D.^2+(tan(theta_i1).*tan(theta_v1).*sin(azimuth_angle)).^2) )./(sec(theta_i1)+sec(theta_v1));
cos_t = max(-1,min(1,cos_t));
clear D h b azimuth_angle;

t = acos(cos_t);
clear cos_t;

O = 1/pi*(t-sin(t).*cos(t)).*(sec(theta_i1)+sec(theta_v1));
O = max(0,O);

k_geo = O - sec(theta_i1)-sec(theta_v1) + 1/2*(1+cos(g_1)).*sec(theta_i1).*sec(theta_v1);
end

