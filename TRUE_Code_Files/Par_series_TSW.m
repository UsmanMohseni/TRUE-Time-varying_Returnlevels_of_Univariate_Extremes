%% Function of Parameter series using Time Sliding window
function [par]=Par_series_TSW(dist_name,Data,MW)

years=size(Data,1);
for i_d=1:years-(MW-1)
    input=Data(i_d:i_d+(MW-1),:);
    nn=input(:);
    nn(isnan(nn))=[];
    
    par_int=eval([dist_name,'fit','(nn)']);
    par(:,i_d)=par_int';
end
end