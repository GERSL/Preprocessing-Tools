function data_nir = TopoCorrectIC(surf_band,index_exclude_cloud,sun_zenith,sun_azimuth, slope,aspect,ndvi)    
    data_nir = surf_band;
    jidim = size(surf_band);
    
    step = 50;
    buffer = 25;
    for i = 1:step:jidim(1)
        for j = 1:step:jidim(2)
            % extract window images
            data_tmp = double(data_nir(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2))));
            data_nir_ori = surf_band(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            index_exclude_cloud_water = index_exclude_cloud(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            sun_zenith_deg = sun_zenith(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            sun_azimuth_deg = sun_azimuth(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            slope_data = slope(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            aspect_data = aspect(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            ndvi_ind = ndvi(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2)));
            ndvi_gt = (ndvi_ind>=0.5&ndvi_ind<=1);
            ndvi_lt = (ndvi_ind<0.5&ndvi_ind>=0);
            % IC topographic correction  
            tmp_ndvi = TopoCorrectICWindow(data_nir_ori,index_exclude_cloud_water,sun_zenith_deg,sun_azimuth_deg, slope_data,aspect_data,ndvi_gt);
            data_tmp(ndvi_gt==1) = double(tmp_ndvi); 
            tmp_ndvi = TopoCorrectICWindow(data_nir_ori,index_exclude_cloud_water,sun_zenith_deg,sun_azimuth_deg, slope_data,aspect_data,ndvi_lt);
            data_tmp(ndvi_lt==1) = double(tmp_ndvi);
            data_nir(max(i-buffer,1):min(i+step+buffer-1,jidim(1)),max(j-buffer,1):min(j+step+buffer-1,jidim(2))) = data_tmp;
        end
    end
end

% CSC+C stratied on cos i with 0.1 increasement. a total of 50,000 pixels.
function [band_data] = TopoCorrectICWindow(band_ori,index_exclude_cloud_water,sun_zenith_deg,sun_azimuth_deg, slope_data,aspect_data,ndvi_ind)
%     History:
%     1. Create this function. (1. January, 2017)
%     2. A total of samples are revised as 40,000 from 50,000. (8. March, 2018)
%     3. When c is calculated as Nan, this function will not make the topo 
%     correction. (8. March, 2018)
   
%     data_nir = data_nir_ori;
%     ori = data_nir_ori;
    [ndvi_ind1] = find(ndvi_ind==true);
    if isempty(ndvi_ind1)
        band_data = [];
        return
    end
    ndvi_ind(1:2:end,1:2:end) = false;
    [ndvi_ind2] = find(ndvi_ind==true);
    loc_ind = ismember(ndvi_ind1,ndvi_ind2)==1;

    
    % extract the corresponding ndvi>0.5 or ndvi<=0.5 pixels
    band_ori = band_ori(ndvi_ind1);
    % make sure the output image size is the same with input
    band_data = band_ori;
    index_exclude_cloud_water = index_exclude_cloud_water(ndvi_ind1);
    sun_zenith_deg = sun_zenith_deg(ndvi_ind1);
    sun_azimuth_deg = sun_azimuth_deg(ndvi_ind1);
    slope_data = slope_data(ndvi_ind1);
    aspect_data = aspect_data(ndvi_ind1);
    
    sun_zenith_rad=deg2rad(double(sun_zenith_deg));
    sun_zenith_cos=cos(sun_zenith_rad);
    sun_zenith_sin=sin(sun_zenith_rad);
    clear sun_zenith_deg sun_zenith_rad sun_zenith_rad;
    
    cos_sita=sun_zenith_cos.*cos(deg2rad(slope_data))+sun_zenith_sin.*sin(deg2rad(slope_data)).*cos(deg2rad(single(sun_azimuth_deg)-aspect_data));
    

    pix_ind = ((double(cos_sita)-double(sun_zenith_cos))<=-0.05)|((double(cos_sita)-double(sun_zenith_cos))>=0.05);
    
    clear aspect_data;
    cos_sita_exclude_cloud=cos_sita(index_exclude_cloud_water&cos_sita>=0&loc_ind&pix_ind);
    data_nir_ori_tmp=band_ori(index_exclude_cloud_water&cos_sita>=0&loc_ind&pix_ind);
    
    if length(cos_sita_exclude_cloud)<=10
        band_data = band_ori;
    else
        c_fitted_1=polyfit(double(cos_sita_exclude_cloud),double(data_nir_ori_tmp),1);
        band_data(pix_ind)=double(band_ori(pix_ind))-c_fitted_1(1,1)*(cos_sita(pix_ind)-sun_zenith_cos(pix_ind));    
    end
end

