function band_corrected = TopoCorrectSCS(band_ori,sun_zenith_deg,sun_azimuth_deg, slope_data,aspect_data)
%TOPOCORRECT makes a topographic correction.
%   Detailed explanation goes here
    band_corrected = band_ori;
    sun_zenith_rad=deg2rad(double(sun_zenith_deg));
    sun_zenith_cos=cos(sun_zenith_rad);
    sun_zenith_sin=sin(sun_zenith_rad);
    clear sun_zenith_deg sun_zenith_rad sun_zenith_rad;
    cos_sita=sun_zenith_cos.*cos(deg2rad(slope_data))+sun_zenith_sin.*sin(deg2rad(slope_data)).*cos(deg2rad(single(sun_azimuth_deg)-aspect_data));
    clear aspect_data;
   
    % when -0.05 < cos_sita - sun_zenith_cos < 0.05, DO NOT make
    % correction. Ref. Tan et al. RSE (2013)
    cor_mask = ((double(cos_sita)-double(sun_zenith_cos))<=-0.05)|...
        ((double(cos_sita)-double(sun_zenith_cos))>=0.05);
    
    band_corrected(cor_mask)=double(band_ori(cor_mask)).*(cos(deg2rad(slope_data(cor_mask))).*sun_zenith_cos(cor_mask))./(cos_sita(cor_mask));
end

