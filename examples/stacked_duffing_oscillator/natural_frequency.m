clear
% close all
NUM_POINTS = 100;
stiff_ratio_21 = linspace(0.1,1,NUM_POINTS);
stiff_ratio_31 = 50;

k1 = 1;
m1 = 10;
m2 = 0.5;
m3 = 1;

% syms k2 k3
freq_ratio = zeros(2,NUM_POINTS);
for iPoint = 1:NUM_POINTS
    k2 = stiff_ratio_21(iPoint)*k1;
    k3 = stiff_ratio_31*k1;
    
    K = [k1 + k2, -k2, 0;
        -k2, k2+k3, -k3;
        0, -k3, k3];
    
    M = eye(3).*[m1;m2;m3];

    [evec,eval] = eig(K,M);

    freq = sqrt(diag(eval));

    freq_ratio(:,iPoint) = freq(2:3)/freq(1);
end

figure;
plot(stiff_ratio_21,freq_ratio)


% NUM_POINTS = 100;
% stiffness_ratio_21 = linspace(0,1,NUM_POINTS);
% stiffness_ratio_31 = linspace(0,1,NUM_POINTS);
% 
% k1 = 1;
% % syms k2 k3
% freq_ratio_21 = zeros(NUM_POINTS,NUM_POINTS);
% freq_ratio_31 = zeros(NUM_POINTS,NUM_POINTS);
% for iPoint = 1:NUM_POINTS
%     for jPoint = 1:NUM_POINTS
%         k2 = stiffness_ratio_21(iPoint)*k1;
%         k3 = stiffness_ratio_31(jPoint)*k1;
% 
%         K = [k1 + k2, -k2, 0;
%             -k2, k2+k3, -k3;
%             0, -k3, k3];
% 
% 
% 
%         [evec,eval] = eig(K);
% 
%         freq = sqrt(diag(eval));
% 
%         freq_ratio_21(iPoint,jPoint) = freq(2)/freq(1);
%         freq_ratio_31(iPoint,jPoint) = freq(3)/freq(1);
%     end
% end
% 
% 
% [X,Y] = meshgrid(stiffness_ratio_21,stiffness_ratio_31);
% figure;
% tiledlayout("flow")
% nexttile
% mesh(X,Y,freq_ratio_21)
% nexttile
% mesh(X,Y,freq_ratio_31)