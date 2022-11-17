%All used PET/SPECT data are publicaly available and downloaded from 
%https://github.com/juryxy/JuSpace/tree/JuSpace_v1.3/JuSpace_v1.3/PETatlas

%Data path
path = '';
%PET/SPECT name
ls = dir(fullfile(path,'data'));
dn = cell(length(ls)-2,1);
for i = 1:length(ls)-2
    dn{i} = ls(i+2).name;
end

%HCP template
[hcp_lh,~,~,HeDEPer1]=y_ReadAll(fullfile(path,'Reslice_lh.HCP.ch2.nii'));
[hcp_rh,~,~,HeDEPer2]=y_ReadAll(fullfile(path,'Reslice_rh.HCP.ch2.nii'));
U_r = unique(hcp_rh(hcp_rh~=0));
U_l = unique(hcp_lh(hcp_lh~=0));

PET_HCP = zeros(size([U_r;U_l],1),length(dn));
PET_HCP_miss = zeros(size([U_r;U_l],1),length(dn));

for jj = 1:length(dn)
    %read PET data
    [df,~,~,~] = y_ReadAll(fullfile(fd,'data',dn{jj}));
    for i = 1:size([U_r;U_l],1)
        if i < 181
            tp_r = df(hcp_rh==i);
            v1 = find(tp_r<=0);
            if isempty(v1) == 0
                tp_r(tp_r<=0) = nan;
                PET_HCP(i,jj) = nanmean(tp_r);
                PET_HCP_miss(i,jj) = length(v1)/length(tp_r);
            else
                PET_HCP(i,jj) = mean(tp_r);
            end
        elseif i > 180
            tp_l = df(hcp_lh==i);
            v2 = find(tp_l<=0);
            if isempty(v2) == 0
                tp_l(tp_l<=0) = nan;
                PET_HCP(i,jj) = nanmean(tp_l);
                PET_HCP_miss(i,jj) = length(v2)/length(tp_l);
            else
                PET_HCP(i,jj) = mean(tp_l);
            end
        end
    end
end

save PETdata_HCP PET_HCP PET_HCP_miss dn;