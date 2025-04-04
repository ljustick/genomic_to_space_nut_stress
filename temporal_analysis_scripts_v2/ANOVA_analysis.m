clear
close all

load DataLL_theta.mat;
load time.mat;

data = DATAlowLat3;
N3 = size(data,3);
data2 = reshape(data,29520,N3);
C = sum(data2,2);

ind = ~isnan(C);%This index is used when back transforming
data3 = data2(ind,:)';
%% 

N = size(data3,2);
SSQ = zeros(N,4);
Coefs = zeros(N,24);
%i=1;
for i = 1:N

    [p,tbl,stats,terms] = anovan(data3(:,i),time_vec(:,1:2),'display','off');

    SSQ(i,:)   = cell2mat(tbl(2:5,2));
    Coefs(i,:) = stats.coeffs;
end

save('ANOVA_output.mat',"Coefs","SSQ","ind")