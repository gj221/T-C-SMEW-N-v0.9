function [data_s_day, data_m_day] = Cflux_post_processing(data, maize_year, soybean_year)

    data_GPP = data.RA_L(2:end,:) + data.NPP_L(2:end,:);


    data_NEE = -data.NPP_L(2:end,:)+(data.R_litter(2:end,:) + data.R_microbe(2:end,:) + data.R_ew(2:end,:));
    
    data_Re = data.R_litter(2:end,:) + data.R_microbe(2:end,:) + data.R_ew(2:end,:)+ data.RA_L(2:end,:);
    
    data_Raut = data.Rmr_L(2:end,:);
    data_Rhet = data.R_litter(2:end,:) + data.R_microbe(2:end,:) + data.R_ew(2:end,:);
    
    data_datum = datevec(data.Date(1:24:end));
    
    % years when maize and soybean in place
    maize_index = ismember(data_datum(:,1),maize_year);
    soybean_index = ismember(data_datum(:,1),soybean_year);
    
    % seperate sim for maize and soybean
    data_GPP_m = data_GPP(maize_index,1);
    data_GPP_s = data_GPP(soybean_index,2);
    
    data_NEE_m = data_NEE(maize_index,1);
    data_NEE_s = data_NEE(soybean_index,2);
    
    data_Re_m = data_Re(maize_index,1);
    data_Re_s = data_Re(soybean_index,2);
    
    data_Rhet_m = data_Rhet(maize_index);
    data_Rhet_s = data_Rhet(soybean_index);
    
    data_Raut_m = data_Raut(maize_index,1);
    data_Raut_s = data_Raut(soybean_index,2);
    
    % annual average daily: [gC/m2/day]
    data_GPP_m_day = zeros(366,1);
    data_GPP_s_day = zeros(366,1);
    
    data_NEE_m_day = zeros(366,1);
    data_NEE_s_day = zeros(366,1);
    
    data_Re_m_day = zeros(366,1);
    data_Re_s_day = zeros(366,1);
    
    data_Rhet_m_day = zeros(366,1);
    data_Rhet_s_day = zeros(366,1);
    data_Raut_m_day = zeros(366,1);
    data_Raut_s_day = zeros(366,1);
    
    j = 1;
    for mm = 1:12
        for dd = 1:31
            
            % maize
            data_m = data_datum(maize_index,:);
            find_day1 = (data_m(:,2)==mm & data_m(:,3)==dd);
            data_GPP_m_day(j) = mean(data_GPP_m(find_day1));
            data_NEE_m_day(j) = mean(data_NEE_m(find_day1));
            data_Re_m_day(j) = mean(data_Re_m(find_day1));
    
            data_Rhet_m_day(j) = mean(data_Rhet_m(find_day1));
            data_Raut_m_day(j) = mean(data_Raut_m(find_day1));
    
            % soybean
            data_s = data_datum(soybean_index,:);
            find_day2 = (data_s(:,2)==mm & data_s(:,3)==dd);
            data_GPP_s_day(j) = mean(data_GPP_s(find_day2));
            data_NEE_s_day(j) = mean(data_NEE_s(find_day2));
            data_Re_s_day(j) = mean(data_Re_s(find_day2));
    
            data_Rhet_s_day(j) = mean(data_Rhet_s(find_day2));
            data_Raut_s_day(j) = mean(data_Raut_s(find_day2));
            
            if sum(find_day1)>0 || sum(find_day2)>0
                j = j+1;
            end
    
        end
    end
    
    data_GPP_s_day(isnan(data_GPP_s_day))=[];
    data_NEE_s_day(isnan(data_NEE_s_day))=[];
    data_Re_s_day(isnan(data_Re_s_day))=[];
    
    data_Raut_s_day(isnan(data_Raut_s_day))=[];
    data_Rhet_s_day(isnan(data_Rhet_s_day))=[];
    
    data_GPP_m_day(isnan(data_GPP_m_day))=[];
    data_NEE_m_day(isnan(data_NEE_m_day))=[];
    data_Re_m_day(isnan(data_Re_m_day))=[];
    
    data_Raut_m_day(isnan(data_Raut_m_day))=[];
    data_Rhet_m_day(isnan(data_Rhet_m_day))=[];

    data_s_day.data_GPP_s_day = data_GPP_s_day;
    data_s_day.data_NEE_s_day = data_NEE_s_day;
    data_s_day.data_Re_s_day = data_Re_s_day;
    data_s_day.data_Raut_s_day = data_Raut_s_day;
    data_s_day.data_Rhet_s_day = data_Rhet_s_day;

    data_m_day.data_GPP_m_day = data_GPP_m_day;
    data_m_day.data_NEE_m_day = data_NEE_m_day;
    data_m_day.data_Re_m_day = data_Re_m_day;
    data_m_day.data_Raut_m_day = data_Raut_m_day;
    data_m_day.data_Rhet_m_day = data_Rhet_m_day;


end