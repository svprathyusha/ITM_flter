function [ output ,ite2] = ITM_new( x)

xo=x;

n=length(x); %length of the image
%% for stopping criteria conditions
%section 2- (C) stpping criteria
e1=1;              % equation (34)    
e2= 2*n^0.5;       %equation (35)
e3=(n-n^0.5)/2;    %equation(36)
e4=n^0.5;          %euation (38)
s4= 0; % for using in the equation (37)
ite_times=0;
ite2=0;
 



while(1)
    g=s4;
    ite2=ite2+1;
    ite_times= ite_times+1;
    %% Outline of the ITM algorithm:
    %step 1 : computing arithmetic mean
    u=mean(x);   %equation(2)
    % step 2 : computing threshold and truncating the input data
  
    % computing threshold first
    % section 2 - (B) Finding dynamic trucation threshold
    % finding higher values than the mean
    xh=x(x>u);  % equation(6)
    nh=length(xh); % size of the higher values than mean
    uh=mean(xh); % mean of the higher values
    dh= uh-u; % computing the higher mean
    %finding lower values than the mean
    xl=x(x<=u); %equation(7)
    nl=length(xl); % size of the lower values than mean
    ul= mean(xl);% mean of the lower values
    dl=u-ul; %computing difference of means
    % first threshold (t1)
    %t = 0.5*(dh+dl); %equation (10) 
   % t= std(x-u);  %second threshold( t2) %equation(11)
    t= mean(abs(x-u)); %third threshold( t3) % euqtion(12)
    %truncating the values 
    bh= u+t; % equation(3)
    bl=u-t;   % equation(3)
    xh_truncated=(x>bh); %finding higher values than mean and threshold
       
    nth=length(xh_truncated); %length of higher truncated values
    x(xh_truncated)=bh; % replacing the higher values with u+t
    xl_truncated=(x<bl); %finding lower values than mean and threshold
     ntl=length(xl_truncated); % length of lower truncated values
    x(xl_truncated)=bl;  % replacing the lower values with u+t
 
   
    %% Section 2 -(C)Stopping criteria
   s1=abs(nh-nl); %equation 34
   s3= abs(nth-ntl); % equation 36
  % s4= abs(nth-ntl); %equation 37
   % equation 38
  if s1<=e1 || ite_times>=e2 || s3>=e3 || ((s3>=e4)&&(s4 == g))% g
   %if s1<=e1 || ite_times>=e2 || s3>=e3 || ((s3>=e4)&&(s3 == s4))
        break;
   end
   s4=s3; %equation 37
end 


   xr=x>bl&x<bh;
        if sum(xr)>(n/2) % to avoid unreliable mean
           output= mean(x(xr)); %equation (5)
        elseif sum(xr)>(n/4)
            output=(mean(x(xr))+mean(x))/2;
        else
            output= mean(x); %equation(5)
        end
 

      