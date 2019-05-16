function [sr_cor,c] = TopoCorrectSCSplusC(sr_ori,clr_mask,sun_zenith_deg,sun_azimuth_deg, slope_data,aspect_data)
                    
% TOPOCORRECTCSCPLUSC SCS + C stratied on cos i with 0.1 increasement.
% a total of 40,000 pixels.
%
% Ref. Qiu, Shi, et al. "Improving Fmask cloud and cloud shadow detection 
% in mountainous area for Landsats 48 images." Remote Sensing of 
% Environment 199 (2017): 107-119.
% 
%
% Inputs:
% sr_ori: orginal surface reflectance
%
% clr_mask: clear sky mask; 1 is clear; 0 is non-clear, such as cloud and
% cloud shadow.
%
% sun_zenith_deg: the zenith angle of sun; unit: decimal degree
%
% sun_azimuth_deg: the azimuth angle of sun; unit: decimal degree
%
% slope: slope derived from DEM; unit: decimal degree
%
% aspect: aspect derived from DEM; unit: decimal degree
%
%
%     History:
%     1. Create this function. (1. January, 2017 by Shi Qiu)
%     2. A total of samples are revised as 40,000 from 50,000. (8. March, 2018 by Shi Qiu)
%     3. When c is calculated as Nan, this function will not make the topo 
%     correction. (8. March, 2018 by Shi Qiu)
%     4. when -0.05 < cos_sita - sun_zenith_cos < 0.05, DO NOT make
%     correction. (27. July, 2018 by Yukun Lin)

    sr_cor = sr_ori;
    
    sun_zenith_rad=deg2rad(double(sun_zenith_deg));
    sun_zenith_cos=cos(sun_zenith_rad);
    sun_zenith_sin=sin(sun_zenith_rad);
    clear sun_zenith_deg sun_zenith_rad sun_zenith_rad;
    cos_sita=sun_zenith_cos.*cos(deg2rad(slope_data))+sun_zenith_sin.*sin(deg2rad(slope_data)).*cos(deg2rad(single(sun_azimuth_deg)-aspect_data));
    clear aspect_data;
    cos_sita_exclude_cloud=cos_sita(clr_mask);
    % random stratied sampling
    cos_sita_min=min(cos_sita_exclude_cloud);
    cos_sita_max=max(cos_sita_exclude_cloud);
    
    % when -0.05 < cos_sita - sun_zenith_cos < 0.05, DO NOT make
    % correction. Ref. Tan et al. RSE (2013)
    cor_mask = ((double(cos_sita)-double(sun_zenith_cos))<=-0.05)|...
        ((double(cos_sita)-double(sun_zenith_cos))>=0.05);
    
% %     % do not use the pixels what we do not need to correct
% %     cos_sita_exclude_cloud = cos_sita_exclude_cloud & cor_mask;
    
%     total_sample=50000;
    total_sample=40000; % enough derived from Fmask 4.0 (Qiu at al., 2018)
    cos_sita_interval=0.1;
    samples_ids= stratiedSampleHanlder(cos_sita_exclude_cloud,cos_sita_min,cos_sita_max,total_sample,cos_sita_interval);
    
    cos_sita_samples=cos_sita_exclude_cloud(samples_ids);
    clear cos_sita_exclude_cloud cos_sita_min cos_sita_max total_sample cos_sita_interval;
    
%     sr_ori_tmp=sr_ori(clr_mask&cor_mask);
    data_samples_sr=sr_ori(samples_ids);
    clear sr_ori_tmp;
    c_fitted=polyfit(double(cos_sita_samples),double(data_samples_sr),1);  
%     figure;plot(double(cos_sita_samples),double(data_samples_nir),'r.');
    c=c_fitted(1,2)/c_fitted(1,1);
    clear c_fitted;
    if isnan(c)
        sr_cor=sr_ori;
        clear sr_ori;
    else
        sr_cor(cor_mask)=double(sr_ori(cor_mask)).*(cos(deg2rad(slope_data(cor_mask))).*sun_zenith_cos(cor_mask)+c)./(cos_sita(cor_mask)+c);
        clear sr_ori;
    end
end

function samples_ids=stratiedSampleHanlder(data_dem_clear,dem_b,dem_t,total_sample,ele_strata)
%STRATIEDSAMPLEHANLDER Further select clear sky (land) pixels to normalize
%BT.
%
% Syntax
%
%     samples_ids=
%     stratiedSampleHanlder(data_dem_clear,dem_b,dem_t,dim,total_sample,ele_strata,distance)
%
% Description
%
%     Stratied sampling method is used by Qiu et al., (2017).
%     History:
%     1. Create this function. (1. January, 2017)
%     2. This stratied sampling method sometimes results in no enough
%     samples. (8. March, 2018)
%
% Input arguments
%
%     data_dem_clear     Clear sky pixels' DEM.
%     dem_b              Basic elevation.
%     dem_t              Highest elevation.
%     dim                Dim for data.
%     total_sample       All clear sky samples (40,000).
%     ele_strata         Stratied sampling along elevation (300 meters).
%     distance           Minmum distance among different samples(450 meters).
%                        Distance rule will be removed if distance = 0.
%     resolution         Spatial resolution (Landsat 30 meters; Sentinel-2 20 meters).
% Output arguments
%
%     Tempd          Nomalized Temperature (BT).
%
%        
% Author:  Shi Qiu (shi.qiu@ttu.edu)
% Date: 8. March, 2018


% %     num_clear_sky_pixels=numel(data_dem_clear);
% %     % no enough clear-sky pixels afther removing the pixels out of the upper
% %     % and lower levels.more than 1/4 total samples (10,000) can be used for 
% %     % estimating lapse rate and c in topo correction.
% %     if num_clear_sky_pixels<total_sample/4 
% %         samples_ids=1:num_clear_sky_pixels;
% %         return;
% %     end

    % compute the number of available strata
    strata_avail=[];
    for i_dem=dem_b:ele_strata:dem_t
        dem_clear_index_tmp=data_dem_clear>=i_dem&data_dem_clear<i_dem+ele_strata;
        if sum(dem_clear_index_tmp)>0
%                 num_strata=num_strata+1;
            strata_avail=[strata_avail;1];
        else
            strata_avail=[strata_avail;0];
        end
        clear dem_clear_index_tmp;
    end
    % equal samples in each stratum
    num_sam_each_strata=round(total_sample/sum(strata_avail));
    clear total_sample;
    samples_ids=[];% to restore samples selected.
    % loop each strata and select samples
    loop_i=1;
    for i_dem=dem_b:ele_strata:dem_t
        if strata_avail(loop_i)
            % find all clear-sky pixels in subgroup.
            samples_ids_tmp=find(data_dem_clear>=i_dem&data_dem_clear<i_dem+ele_strata);
            % randomly selection.
            samples_ids_rand=samples_ids_tmp(randperm(numel(samples_ids_tmp))); 
            clear samples_ids_tmp;
            num_tmp=size(samples_ids_rand,1);
            num_max_tmp=num_sam_each_strata;
            if num_tmp>num_max_tmp
                num_tmp=num_max_tmp;
            end
            clear num_max_tmp;
            samples_ids_rand_tmp=samples_ids_rand(1:num_tmp);
            clear num_tmp;
            % store data
            samples_ids=[samples_ids;samples_ids_rand_tmp];
            clear samples_ids_rand samples_ids_rand_tmp;
        end
        loop_i=loop_i+1;
    end
    clear samp_dis num_sam_each_strata;

% %     % no enough clear-sky pixels afther removing the pixels out of the upper
% %     % and lower levels.more than 1/4 total samples (10,000) can be used for 
% %     % estimating lapse rate and c in topo correction.
% %     if numel(samples_ids)<total_sample/4
% %         samples_ids=1:num_clear_sky_pixels;
% %         return;
% %     end
end