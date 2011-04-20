function [ilsTime] = frame(sat_number, epoches, newILS)

%function [] = raim1acsd1(p_runs,test) 
% DESCRIPTION: RAIM with the sd accumulated linear model. 
%              p_runs - number of runs 
%              test   - test statistic (optional parameter): 
%                 1 - Maximum Residual test statistic, 
%                 2 - Sturza test statistic, 
%                 3 - Brenner test statistic, 
%                 4 - all test statistics (default value); 
 
% ------------------------------------------------------------------------------- 
% Constants 
% ------------------------------------------------------------------------------- 
DIM		= 3; 
RAD_P_DEG 	= pi/180; 	% radians per degree and vice-versa; 
DEG_P_RAD 	= 1 / RAD_P_DEG; 
L1_F 		= 1575.42e6; 	% L1 carrier frequency 
C 		= 3e+8; 	% speed of light 
L1_W 		= C/L1_F; 	% L1 wavelength 
EARTH_RADIUS 	= 6378137; 	% earth radius 

% ------------------------------------------------------------------------------- 
% Parameters 
% -------------------------------------------------------------------------------
intflag = 0;
%%%%
global CC_m;
%%%%

dout=0; %display plots (the last p_runs) 
pout=0; %print plots (the last p_runs) 
%t_init=rand(1)*5000+26000; 
t_init =30000;     %initial time (offset for epoch),26000<= t_init <=31000,from Hunzinger thesis  
epoch_max= epoches;     %maximum run time  
noise_sigma=0.001; 
p_base=100; 
p_unc=1000; 
p_speed=100;
min_elev 	= 5*RAD_P_DEG; % 5 degrees, minimum elevation for visibility 
t_epoch 	= 1;		   % in seconds 
n_vis_sat_max	= sat_number;     % must be:  n_vis_sat_max >= 7 
ilsTime = zeros(1,epoches);
n_raim=7;                  % must be:  5 <= n_raim <= n_vis_sat_max 
k_epochs=4;                % the number of accumulated epochs for the lin. model 
 
% ------------------------------------------------------------------------------ 
p_ephd=yuma2mat;    % load and convert GPS data 
 
 
%------------------------------------------- 
 % satellite related information extracted from ephemeris data 
 [n_sat,ephpars] = size(p_ephd); % no. of satellites, ephemeris parameters for YUMA 
 
 % receiver positions  
 baselength 	= p_base;		% initial baseline length (m) 
 % receiver 1 is static 
 x_rcvr1 	= [EARTH_RADIUS 0 0];	% position of static reciever 1 
 
 % receiver 2 motion model 
 x_rcvr2_r       = baselength*sqrt(2)/2; 
 x_rcvr2_h       = baselength*sqrt(2)/2; 
 s_rcvr2		 = p_speed;	% m/s 
 
 % receiver 2 initial position 
 x_rcvr2 = rcvrmdl(x_rcvr2_h,x_rcvr2_r,s_rcvr2,t_init*t_epoch); 
 
 % initial guess at position and velocity 
 x_unc = [p_unc p_unc p_unc]; 
 
 % position uncertainty (provided by pseudo-range) 
 dx_rcvr2_est 	= (x_rcvr2-x_rcvr1)+x_unc; 
 %%%%%%%%%%%%%%%%%%%%%%%% 
 %CC_alg Initial x1, for here we use the same estimate x as 1'st one(vg.lee) 
 %%%%%%%%%%%%%%%%%%%%%%%% 
 CC_x1  = dx_rcvr2_est'; 
 
 % ------------------------------------------------------------------------------ 
 % Initialization (Observations and Measurements) 
 % ------------------------------------------------------------------------------ 
 % find inital the visible satellites to 2 receivers, the single difference, and the integer ambiguity at the start epoch
 [visible sd iNs] = KFobserve_sd_InterAmbiguity(x_rcvr1,x_rcvr2,dx_rcvr2_est,noise_sigma,t_init*t_epoch,min_elev,p_ephd,n_vis_sat_max); 
 n_vis_sat=length(visible); 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %CC_alg Initial y1 
 %CC_alg Initial Number of satelite 
 %here we use the same initial y and m-number of satelites 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 CC_y1 = sd'; 
 CC_m = n_vis_sat;
  enew                        = ones(CC_m-1,1);
J = [enew,-eye(CC_m-1)];
F = eye(CC_m-1) - (enew*enew')/(CC_m - sqrt(CC_m));
 %-----------receiver clock bias--------------------------- 
 rclbias=zeros(2,1); 
 Rcl=rcfac(t_epoch); 
 % ------------------------------------------------------------------------------ 
 % Initialization (Design Matrix) 
 % ------------------------------------------------------------------------------ 
 % form the design and geometry matrices 
 % 
 %Initial E, 
 %got it trough initial x 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 E = geometry(x_rcvr1,dx_rcvr2_est,visible,p_ephd,t_init*t_epoch); 
 CC_E1 = E/L1_W; 
 
 % ------------------------------------------------------------------------------ 
 % Simulation Main Loop 
 % ------------------------------------------------------------------------------
 epoch=0;
 t_fault=0;
 p_failsat=0;
 p_break=0;
 
%%% compute P 
 e1 = zeros(CC_m,1);
 e1(1,1) = 1;
 e = ones(CC_m,1);
 u = e1 - e / sqrt(CC_m);
 P = eye(CC_m) - u*(2/(u'*u))*u';
 P2 = P(:,2:CC_m);
 
 x_j = CC_x1;
 E_j = CC_E1;
 y_j = CC_y1;
 
 flag=0;
 %%%
 error_vector = zeros(epoch_max,1);
 cov_d_vector = zeros(epoch_max,1);
 cov_x_vector = zeros(epoch_max,1);
while p_break==0, 
   epoch=epoch+1;
   % calculate the actual receiver position 
   % model the roving receiver as constant velocity model 
   x_rcvr2 = rcvrmdl(x_rcvr2_h,x_rcvr2_r,s_rcvr2,epoch*t_epoch);    
   % find the current time visible satellites and single difference at current epoch
   [visible sd] = KFobserve_sd(x_rcvr1,x_rcvr2,dx_rcvr2_est,iNs,noise_sigma,(epoch+t_init)*t_epoch,min_elev,p_ephd,n_vis_sat_max);
   %-----------receiver clock bias--------------------------- 
   rclbias=rclock(rclbias,t_epoch,Rcl); 
   sd=sd-C*rclbias(1)/L1_W;    
   %for each epoch in CC_alg, 
   %we have to get y 
   CC_y = sd'; 
   %---------------------------------------- 
   %    POSITION, YOU HAVE TO START HERE 
   %----------------------------------------  
   
   %get time instant of flight and next y
   t_curr = (t_init+epoch)*t_epoch;
   y_j_1 = CC_y;
   
   % if S_j is not full column rank
   if flag==0   
   
      E_j_1 = E_j;
     
      %computer U_j,W_j, and R_j, assume that P2'*Ej has full column rank 
      [U_j,W_j,R_j] = transform(P2'*E_j);
     
      %compute x_j_1
      x_j_1 = x_j + inv(R_j)* U_j*P2'*(y_j_1 - y_j);

      %use x_j_1 to construct E_j_1;
      E_j_1 = geometry(x_rcvr1,x_j_1',visible,p_ephd,t_curr);
      E_j_1 = E_j_1/L1_W;
      
      %compute U_j,W_j, and R_j, assume that P2'*Ej_j_1 has full column rank 
      [U_j_1,W_j_1,R_j_1] = transform(P2'*E_j_1);
      
      %use U_j_1, R_j_1 to get new proved estimate 
      x_j_1 = inv(R_j_1)*U_j_1*(P2'*E_j*x_j+P2'*(y_j_1-y_j));

      %construct the new E_j_1 using proved new estiamte x_j_1 
      E_j_1 = geometry(x_rcvr1,x_j_1',visible,p_ephd,t_curr);
      E_j_1 = E_j_1/L1_W;
      
   else % else S_j is full column rank
      E_j_1 = geometry(x_rcvr1,x_j',visible,p_ephd,t_curr);
      E_j_1 = E_j_1/L1_W;
   end
   
   %compute U_j_1,W_j_1 and R_j_1 assuming that P2'E_j_1 is full column rank
   [U_j_1,W_j_1,R_j_1] = transform(P2'*E_j_1);
   
   %get w_j_1 and w_j_1
   u_j_1 = U_j_1*(P2'*y_j_1);
   w_j_1 = W_j_1 * (P2'*y_j_1);
      
  
   % call function HQR to compute S_j_1 and w_cap_j_1
   if epoch == 1 % if epoch is 1
       [S_j_1, w_cap_j_1] = HQR(W_j_1*F, w_j_1);
   else % else use prevoius results to get new results
       [S_j_1, w_cap_j_1] = HQR([S_j; W_j_1*F], [w_cap_j; w_j_1]);
   end

   % check if S_j_1 is full colum rank
   if flag==0
      [m,n]=size(S_j_1);
      if rank(S_j_1) == n % if S_j_1 is full column rank, set flag be 1
         flag=1;
         start_epoch=epoch;
      end
   end
   
   % the following code is only for full column rank (that is the part 3.4 & 3.5 in paper)
   if flag == 1
          

        if intflag == 0
         [R,Z,y] = reduction(S_j_1,w_cap_j_1);
               noise = noise_sigma/L1_W;
         sizeR = size(R,2);
            %sizeSj2 = size(S_j_12,2);
            result = 1;
            for i = 1:sizeR     
                entry = R(i,i)/(2*noise);
                result = result*(2*normcdf(entry,0,1) - 1);
            end
           
     
           if(result > 0.9) %1e-8)
               
                disp('FIXED');
                tic;
                INT = ils(S_j_1,w_cap_j_1,1,newILS);
                ilsTime(epoch) = toc;
                d_j_1 = INT;
                intflag = 1;
                fix_epoch = epoch
                % fix_result = res(epoch);
                REAL = S_j_1\w_cap_j_1
                ESTDDIA = d_j_1
                TrueDDIA = [J*iNs']
                
            else
                d_j_1 = S_j_1\w_cap_j_1;
            end
     
       else
          d_j_1 = INT;
       end % end if intflag == 0
   
      % compute convariance only after flag is 1
      %%compute the convariance of x_j_1
      temp = [R_j_1, U_j_1; zeros(CC_m-1,3),S_j_1];
      %% call function Given_transform to get R_bar_j_1 
      R_bar_j_1 = Given_transform(temp);
      cov_x = noise_sigma^2*inv(R_bar_j_1'*R_bar_j_1);
      trace_cov_x = trace(cov_x);
      cov_x_vector(epoch) = trace_cov_x;
       %compute x_j_1 using equation (15)
      x_j_1 = R_j_1\(u_j_1-U_j_1*F*d_j_1);
   end
        
   % compute the position of roaming receiver at the current epoch
   x=x_j_1+x_rcvr1';
   error = norm(x-x_rcvr2');
   %% compute the difference between the true position and estimated position
   error_vector(epoch)=error;
       
   % update variables for next step in while loop   
   E_j = E_j_1;
   x_j = x_j_1;
   y_j = y_j_1;
   S_j = S_j_1;
   w_cap_j = w_cap_j_1;
   
   % test the while loop temination
   if epoch == epoch_max
       p_break = 1;
   end;
 
end; %epochs 


%plot the figures
x = 0:.1:10;
semilogy(x,10.^x,'w');
xlabel('epoch (s)');
ylabel('error in position estimate(m)');
grid;
title(['satellite number is ',int2str(n_vis_sat_max)])
hold on
plot(1,error_vector(1),'r');
for i=2:epoch_max
   x = (i-1):0.001:i;
   y = (error_vector(i) - error_vector(i-1))*(x-i+1)+error_vector(i-1);
   plot(x,y,'r');
end
hold off

figure;
x = 0:.1:10;
semilogy(x,10.^x,'w');
xlabel('epoch (s)');
ylabel('covariance of the error in position estimates (m)');
grid;
title(['satellite number is ',int2str(n_vis_sat_max)])
hold on
plot(start_epoch,cov_x_vector(2),'r');
for i=start_epoch+1:epoch_max
   x = (i-1):0.001:i;
   y = (cov_x_vector(i) - cov_x_vector(i-1))*(x-i+1)+cov_x_vector(i-1);
   plot(x,y,'r');
end
hold off

