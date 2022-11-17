%load PET/SPECT data
load('PETdata_HCP.mat')
%permutation 5000 times
permN = 5000;
%%
%load brain statistic values permuted 5000 times
load('t_perm_5000_area_l5.mat')

%percentage of missing: allowed missing rate in each HCP region
pert = 0.5;
%region size
n = zeros(length(dn),1);
%spearman's correlation
sp_r_area = zeros(length(dn),1);
sp_p_area = zeros(length(dn),1);
sp_p_area_perm = zeros(length(dn),1);

for i = 1:length(dn)
    rd = find(PET_HCP_miss(:,i)<= pert);
    n(i) = length(rd);
    pet = PET_HCP(rd,i);
    %real t value
    t_area = t_area_l5(rd);
    [sp_r_area(i),sp_p_area(i)] = corr(pet,t_area,'type','Spearman');
    
    %permutation
    r_area_p = zeros(permN,1);
    parfor jj = 1:permN
        t_area_p = t_area_l5_p5000(rd,jj);
        [r_area_p(jj),~] = corr(pet,t_area_p,'type','Spearman');
    end
    sp_p_area_perm(i) = length(find(r_area_p>=abs(sp_r_area(i))))/permN;
end

%bootstrapping
bootN = 10000;
sp_r_boot = zeros(length(dn),bootN);

for i = 1:length(dn)
    rd = find(PET_HCP_miss(:,i)<= pert);
    pet = PET_HCP(rd,i);
    t_area = t_area_l5(rd);
    parfor j = 1:bootN
        rd_b = randsample(length(rd),length(rd),true);
        pet_b = pet(rd_b);
        t_b = t_area(rd_b);
        [sp_r_boot(i,j),~] = corr(pet_b,t_b,'type','Spearman');
    end
end