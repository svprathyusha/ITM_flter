clc;
clear;
close all;
%% creating constant image
im=100*ones(100,100 ); % constant image

[m,n]=size(im); % size of image
x=10;% A total of x experiments
N = [3 5 7 9]; % filter window sizes

im_mean=zeros(m,n);
im_median=zeros(m,n);
ITM1=zeros(m,n);
ITM2=zeros(m,n);
ITM3=zeros(m,n);

ite_time_array = zeros(x,4);
ite_new_array = zeros(x,4);

MAE_mean_all=zeros(1,length(N));
MAE_median_all=zeros(1,length(N));
MAE_ITM1_all=zeros(1,length(N)); 
MAE_ITM2_all=zeros(1,length(N));
MAE_ITM3_all=zeros(1,length(N));

for ii=1:x
    
    ii

%% noise
% gausian noise
gau_noise=normrnd(0,15,m,n); gau_noise=gau_noise-mean2(gau_noise);normal_dist= rand(m,n);
sigma=std2(gau_noise);
 
%% gaussian+laplacian noise
% laplacian noise
pd = makedist('laplace','mu',0,'sigma',1.3*(sigma));
lap_noise= random(pd,m,n);

b_mix = 0.50 ;
im_noise= (normal_dist<=b_mix).*gau_noise+(normal_dist>b_mix).*lap_noise+im;



value = mean((im_noise(:)-mean(im_noise(:))).^4)/mean((im_noise(:)-mean(im_noise(:))).^2).^2
%% setting filter window
    for jj=1:4
        N1=(N(jj)-1)/2; % filter window
 %% image filtering 
         for i=1:m
            for j=1:n
                %setting boundaries
                i_arr = max(1,i-N1):min(i+N1,m); % boundaries
                j_arr = max(1,j-N1):min(j+N1,n); % boundaries
                block=im_noise(i_arr,j_arr);
                im_mean(i,j)= mean(block(:)); % mean filter
                im_median(i,j)=median(block(:)); % median filter
          
                [ITM1(i,j),ite_times]=ITM_filter(block(:),1); % iterative truncated mean filter 1
               
                [ITM2(i,j),ite_times]=ITM_filter(block(:),2); % iterative truncated mean filter 2 
                
                [ITM3(i,j),ite_new]=ITM_new(block(:));
                
                ite_time_array(ii,jj)=ite_times; 
                ite_new_array(ii,jj)=ite_new; 
                
            end    
         end
  %% Mean Absolute error 

        MAE_mean(jj)=mean2(abs(im-im_mean));
        MAE_median(jj)=mean2(abs(im-im_median));
        MAE_ITM1(jj)=mean2(abs(im-ITM1));
        MAE_ITM2(jj)=mean2(abs(im-ITM2));
        MAE_ITM3(jj)=mean2(abs(im-ITM3));
         
        MAE_mean_all(jj)=MAE_mean_all(jj)+MAE_mean(jj)/x;
        MAE_median_all(jj)=MAE_median_all(jj)+MAE_median(jj)/x;
        MAE_ITM1_all(jj)=MAE_ITM1_all(jj)+MAE_ITM1(jj)/x;
        MAE_ITM2_all(jj)=MAE_ITM2_all(jj)+MAE_ITM2(jj)/x;
        MAE_ITM3_all(jj)=MAE_ITM3_all(jj)+MAE_ITM3(jj)/x;

    end

end %to end the x loop
%% plotting 
 
%figure,
% plot(N.^2,MAE_mean_all./MAE_median_all,'y-d',N.^2,MAE_ITM1_all./MAE_median_all,'r-o', ...
% N.^2,MAE_ITM2_all./MAE_median_all,'g-s',N.^2,MAE_ITM3_all./MAE_median_all,'k-*',N.^2,MAE_median_all./MAE_median_all,'b-*')
% legend('mean','ITM1','ITM2','ITM3','median')
%  
% set(gca,'XTick',[9 25 49 81] );
% xlabel('filter size n')
% ylabel('MAE normalized by that of median')
% % %ylim([0.99 1.2])
 figure(1),  plot(N.^2,MAE_ITM3_all./MAE_median_all,'k-*', ...
  N.^2,MAE_ITM1_all./MAE_median_all,'r-o',N.^2,MAE_ITM2_all./MAE_median_all,'g-s')
  legend('ITM3','ITM1','ITM2')
 set(gca,'XTick',[9 25 49 81] );
 xlabel('filter size n')
 ylabel('MAE normalized by that of median')
% ITM_filter_avg= mean(ite_time_array)
% ITM_new_avg= mean(ite_new_array)
% figure,
% subplot(1,2,1);imshow(uint8(im_noise));title('noisy image');
% subplot(1,2,2);imshow(uint8(ITM1));title('ITM1 filter');
% figure,
% subplot(1,2,1);imshow(uint8(ITM2));title('ITM2 filter');
% subplot(1,2,2);imshow(uint8(ITM3));title('ITM3 filter');

