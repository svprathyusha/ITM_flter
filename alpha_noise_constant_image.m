clc;clear;close all;
im=100*ones(100,100);

[m,n]=size(im);
 N=5;
 x=1;
 N1=(N-1)/2;
 ite_time_array = zeros(x,1);
  ite_new_array = zeros(x,1);
  
 im_mean=zeros(m,n);
 im_median=zeros(m,n);
 ITM1=zeros(m,n);
 ITM2=zeros(m,n);
 ITM3=zeros(m,n);

 MSE_mean_all=zeros(1,length(N));
MSE_median_all=zeros(1,length(N));
MSE_ITM1_all=zeros(1,length(N)); 
MSE_ITM2_all=zeros(1,length(N));
MSE_ITM3_all=zeros(1,length(N));


for ii=1:x
 ii
%% noise
 %alpha stable noise
 noise = makedist('Stable','alpha',0.5,'beta',0,'gam',10,'delta',0);
alpha_noise=random(noise,m,n);
 im_noise=alpha_noise+im;
 %cliping
  im_noise(im_noise>255)=255;
  im_noise(im_noise<0)=0;
  
value = mean((im_noise(:)-mean(im_noise(:))).^4)/mean((im_noise(:)-mean(im_noise(:))).^2).^2
 for i=1:m
    for j=1:n
       
        %setting boundaries
        i_arr = max(1,i-N1):min(i+N1,m);
        j_arr = max(1,j-N1):min(j+N1,n);
         block=im_noise(i_arr,j_arr);
        im_mean(i,j)= mean(block(:)); % mean filtering
        im_median(i,j)=median(block(:)); % median filtering 
      
        [ITM1(i,j),ite_times]=ITM_filter(block(:),1); %ITM1 filter
        [ITM2(i,j),ite_times]=ITM_filter(block(:),2); % ITM2 filter
        [ITM3(i,j),ite_new]=ITM_new(block(:)); % ITM2 filter
      
       ite_time_array(ii)= ite_times;
       ite_new_array(ii)= ite_new;
       

    end
 end

MSE_mean=immse(im,im_mean);
MSE_median= immse(im,im_median);
MSE_ITM1= immse(im,ITM1);
MSE_ITM2= immse(im,ITM2);
MSE_ITM3= immse(im,ITM3);

 MSE_mean_all=MSE_mean_all+MSE_mean/x;
 MSE_median_all=MSE_median_all+MSE_median/x;
 MSE_ITM1_all=MSE_ITM1_all+MSE_ITM1/x;
 MSE_ITM2_all=MSE_ITM2_all+MSE_ITM2/x;
 MSE_ITM3_all=MSE_ITM3_all+MSE_ITM3/x;
 
 end   
 MSE_mean_all
 MSE_median_all
 MSE_ITM1_all
 MSE_ITM2_all
  MSE_ITM3_all
 
avg_ite=mean(ite_time_array)
avg_ite_new=mean(ite_new_array)

% figure, imshow(uint8(im))
% figure, imshow(uint8(im_noise))
% figure, imshow(uint8(im_mean))
% figure, imshow(uint8(im_median))
% figure, imshow(uint8(ITM1))
% figure,imshow(uint8(ITM2))
