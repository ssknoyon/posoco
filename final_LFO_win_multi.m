 clc
clear all
format short

display('Enter samples in new6.xls format');
%input('\n Press enter when done ');

w1 = 2*pi*50;
% t = input('\n Enter the sampling time (std PMU Ts is 1/25 ) = ');     %1/25;                   % Sampling time/
t=1/25;

y_in = xlsread('voltage(agartala 6th sept).xlsx','B100:B1500');     % input file name 'new1' in xls format of Dadri data

y = y_in/1000000;

W = size(y,1);           % size of total samples (min=200, max = 800-1000)


if (W<200)
    display('No of Samples must be atleast 200; Enter samples again to new6.xls');
    return
end

Cit = size(y,2);

N = 500;%input('\n Enter the window size (if not known type 0; default is 100) = '); %N = 40;                   
                          % window size in consideration 
                          % min length = length of the lowest freq needed
  %N = W/5;                       
if (N < 100 || mod(N,10)~= 0)
    N = 100;
end

% -------------------------------------------------------------------------
% WARNING: choose window size such tat minimum freq is also involved 
%          (atleast 200 samples per window) 
% -------------------------------------------------------------------------

% p = input('\n Enter the order or no of roots desired (if unknown type 0)= ');
p = 100;

if (p == 0)
    p = floor(N/2);
end

% p = 20;                     % order = p, 
                             %take care that order is small preferably<=10                     



 p1 = N;             % if initial value of next window is the value after the final value of prev window
%  p1 = p;            %if window is based on moving window algorithm

% perc = input('\n Enter the percent reduction of singular matrix (if not known type 0) = ');
perc = 1;

if (perc == 0)
    perc = 1;
end

Y1 = zeros(1,N-p);
Y = zeros(1,N);
Y2 = zeros(N-p,p);
B = zeros(N,p);
zr = zeros(1,p);
f = zeros(1,p);
d = zeros(1,p);
sig = zeros(1,p);
w = zeros(1,p);
AMP = zeros(1,p);
Ph = zeros(1,p);


for i1 = 1:Cit
    Q = 0;
    Q1 = 0;
    k = 0;                         % window no initialisation
    r = 0;                       % excel file row initialisation
    x1 = 0; 
for i = 1:W
  s = p+1+k;
  if(s<=W-N+p+1)  
    for j=1:(N-p)
       Y1(j) = y(s,i1);                 % order (N-p)x1
       s=s+1;
    end
    
    for j=1:N
        Y(j) = y(k+j,i1);
    end
    
    x = p+1+k;
    for j=1:(N-p)
       for m=1:p
          Y2(j,m) = y(x-m,i1);           % order (N-p)x(p)
       end
    x=x+1;
    end
    
    %-----------------------
    
    [C,S,D]=svd(Y2);                  % SVD of Y2 matrix 
                                      % it gives S in increasing order of
                                      % eigen values
      
                                                                       
    [S_R,S_C]=size(S);                % S_R = no of rows of S
                                      % S_C = no of cols of S
    if(S_R>S_C)
        for j=1:(S_R - S_C)
            S = S(1:S_C,1:S_C);      % if no of rows is greater than cols 
                                     % delete the rows by making the matrix 
                                     % a square matrix of size S_C
                                     
            C = C(1:S_R,1:S_C);      % delete the corresponding cols in C
        end
    end 
    if(S_C>S_R)                      % same as above when cols>rows  
        for j=1:(S_C - S_R)          % just tat corr rows in D' is deleted 
            S = S(1:S_R,1:S_R);
            D = D(1:S_C,1:S_R);
        end
    end
    
    
    
    red_size = 0;                        % initialise reduced order
    for j=1:p                            
        if(S(j,j)>=perc*S(2,2)/100)
            red_size = red_size + 1;    % if the eigen value is less than  
                                        % 1% of max eigen value which is 
                                        % the first one in S then stop
                                        % updating the counter red_size
        end
    end
    
        S = S(1:red_size,1:red_size);   % make S to the order red_size
        C = C(1:S_R,1:red_size);        % delete corr cols in C  
        D3 = D';
        D2 = D3(1:red_size,1:S_C);      % delete corr rows in D'
        D = D2';
            
    
    Y3 = D*inv(S)*C'*Y1';    % Y3=pinv(Y2)*Y1' can also be written this way 


    
    c1 = [1;-Y3];
    z = roots(c1);
    
    
    %-----------------------
    

    Z = z';
   

    for j=1:N
       for m=1:p
          B(j,m)=z(m)^(j-1);
       end
    end
    
    %--------------------
    
    [C1,S1,D1]=svd(B);          % do the same way as done before for pinv
    [S_R1,S_C1]=size(S1);
    if(S_R1>S_C1)
        for j=1:(S_R1 - S_C1)
            S1 = S1(1:S_C1,1:S_C1);
            C1 = C1(1:S_R1,1:S_C1);
        end
    end 
    if(S_C1>S_R1)
        for j=1:(S_C1 - S_R1)
            S1 = S(1:S_R1,1:S_R1);
            D1 = D1(1:S_C1,1:S_R1);
        end
    end
    
    
    u = D1*inv(S1)*C1'*Y';
    
    %--------------------


    for j=1:p
        zr(j)=log(Z(j));
        f(j)=imag(zr(j))/(2*pi*t);
        d(j)=real(zr(j))/t;
       sig(j) = d(j);
       w(j) = f(j);
       AMP(j) = 2*abs(u(j));
       Ph(j) = angle(u(j));
    end
    
    x = 1;
    x2 = 1;
    

    for j=r+1:r+p
            Q(j,1) = x1+1;
            Q(j,2) = sig(x);
            Q(j,3) = w(x);
            Q(j,4) = AMP(x);
            Q(j,5) = Ph(x)*180/pi;
            Q(j,6) = AMP(x)*exp(sig(x)/w(x));
            x=x+1;
            
            M(j,1) = x1+1;
            M(j,2) = Z(x2);
            M(j,3) = u(x2);
            
            x2 = x2+1;
    end

    r=r+p;     % change this 

    k = k + p1;        % also change this ------------- (1)
    
    x1 = x1+1;
  end

end



% xlswrite('out_win',Q,i1);      % output in xls format with name out_win


% ----------------------------------------------------------------------
% MODE CLASSIFICATION
% ----------------------------------------------------------------------

r=0;
h=0;
o=0;

%f1 = input('\n Enter the largest frequency value to be output (press 0 for default = 2.0) =  '); 
f1 = 2;
if (f1==0)
    f1 = 2.0;
end

for i=1:Q(size(Q,1),1)
    j=1;
    
    for k=r+1:r+p
        if(Q(k,3)>0  && Q(k,3)<f1 && Q(k,4)>0.000001*max(Q(r+1:r+p,4))) 
            R(j,1) = Q(k,2);
            R(j,2) = Q(k,3);
            R(j,3) = Q(k,4);
            R(j,4) = Q(k,5);
            R(j,5) = Q(k,6);
            j=j+1;
        end
    end
    
    R = sortrows(R,-5);
    r1 = size(R,1);
    Q1(h+1:h+r1,1:5) = R(1:r1,1:5);
    Q1(h+1:h+r1,6) = o+1;
    o = o+1;
    h=h+r1;
    r=r+p;
end
R2 = Q1;
n_o_w = x1;
R_tot = size(R2,2);

Q101 = {'Damping','Frequency','Amp','Phase','Factor','Window no'};

xlswrite('sorted_single',Q101,i1,'A1');
xlswrite('sorted_single',Q1,i1,'A2');   


% subplot(2,2,3) implies two sections vertically, horizontally and the 
% 3rd section is in consideration 



ox = floor(W/N);

ox1=0;


% for j=1:ox
%    j2=ox1+1:ox1+N;
%    subplot(1,ox,j);
%    axis([ox1+1 ox1+N ymin ymax]);
%    plot(j2, y(j2),'b',j2,o(j2),'r');
%    ox1 = ox1 + p1;     % change when (1) changes
% end
x4=0;
for j=1:x1
    x2=1;
    x3=0;
    
    for k=1:(size(Q1,1))
        if(Q1(k,6)==x2)
            x3 = x3+1;
        end
    end
    
    
    
    for j3=1+ox1:N+ox1
        for k =1+x4:x3+x4
            o3(j3)=0;
        end
    end
    
    
    for j3=1+ox1:N+ox1
        for k =1+x4:x3+x4
            o3(j3) = o3(j3)+(0.5*Q1(k,3)*exp(Q1(k,1)*0.04*(j3-1-ox1))*cos((2*pi*0.04*(j3-1-ox1)*Q1(k,2)) + Q1(k,4)));
        end
    end
    
    
    ox1 = ox1 + p1; 
    
    x2 = x2+1;
    x4 = x4+x3;
end

ox1=0; 
for j=1:x1
%    figure(i1)
%    subplot(2,x1,j) 
   j2=1+ox1:N+ox1;
%    plot(j2, y(j2,i1),'b');
%    xlabel('Sample number');
%    ylabel('Voltage');
%    grid on
%    subplot(2,x1,j+x1)
%    plot(j2,o3(j2),'r');
%    xlabel('Sample number');
%    ylabel('Mode amplitude');
%    grid on
   ox1=ox1+p1;
end

end

display('------------------------------------------------------------------');
display('Check the excel files sorted_single and sorted_multi for modes and sorting');
display('------------------------------------------------------------------');

% ------------------------------------------------------------------------
% If overlap of figures are necessary then include zero frequency mode
% but see to it that no more than 1 mode with zero freqeuency mode is
% present 
% ------------------------------------------------------------------------


%input('\n Press Enter to continue to multi signal analysis'); 

% ------------------------------------------------------------------------
% MULTI SIGNAL ANALYSIS
% ------------------------------------------------------------------------

k = 0;               % window no initialisation
r = 0;               % excel file row initialisation
x1 = 0;              % roots and residues window no counter initialisation
k1 = 1;

for i = 1:W
  s = p+1+k;                 % value to start from in Y1
  if(s<=W-N+p+1)  
    for j=1:(N-p)*Cit           
       Y1(j) = y(s,k1);      % order (N-p).cx1
       k1 = k1+1;
       if (k1>Cit)
           s=s+1;
           k1 = 1;
       end
    end
    
    for j=1:N
        for l=1:Cit
            Y(j,l) = y(k+j,l);
        end
    end
    
    x = p+1+k;   % value to start from in Y2
    k1 = 1;
    for j=1:(N-p)*Cit
       for m=1:p
          Y2(j,m) = y(x-m,k1);           % order (N-p)x(p)
       end
       k1 = k1+1;
       if (k1>Cit)
           x = x+1;
           k1 = 1;
       end
    end
    
    %-----------------------
    
    [C,S,D]=svd(Y2);                  % SVD of Y2 matrix 
                                      % it gives S in increasing order of
                                      % eigen values
      
                                                                       
    [S_R,S_C]=size(S);                % S_R = no of rows of S
                                      % S_C = no of cols of S
    if(S_R>S_C)
        for j=1:(S_R - S_C)
            S = S(1:S_C,1:S_C);      % if no of rows is greater than cols 
                                     % delete the rows by making the matrix 
                                     % a square matrix of size S_C
                                     
            C = C(1:S_R,1:S_C);      % delete the corresponding cols in C
        end
    end 
    if(S_C>S_R)                      % same as above when cols>rows  
        for j=1:(S_C - S_R)          % just tat corr rows in D' is deleted 
            S = S(1:S_R,1:S_R);
            D = D(1:S_C,1:S_R);
        end
    end
    
    
    
    red_size = 0;                        % initialise reduced order
    for j=1:p                            
        if(S(j,j)>=perc*S(2,2)/100)
            red_size = red_size + 1;    % if the eigen value is less than  
                                        % 1% of max eigen value which is 
                                        % the first one in S then stop
                                        % updating the counter red_size
        end
    end
    
        S = S(1:red_size,1:red_size);   % make S to the order red_size
        C = C(1:S_R,1:red_size);        % delete corr cols in C  
        D3 = D';
        D2 = D3(1:red_size,1:S_C);      % delete corr rows in D'
        D = D2';
            
    
    Y3 = D*inv(S)*C'*Y1';    % Y3=pinv(Y2)*Y1' can also be written this way 


    
    c1 = [1;-Y3];         % coefficients of the eqn is negative of Y3
    z = roots(c1);        % roots of the equation gives freq and damping
    
    
    %-----------------------
    

    Z = z';
   

    for j=1:N
       for m=1:p
          B(j,m)=z(m)^(j-1);   % matrix of roots powers, order = (Nxp)
       end
    end
    
    %--------------------
    
    [C1,S1,D1]=svd(B);          % do the same way as done before for pinv
    [S_R1,S_C1]=size(S1);
    if(S_R1>S_C1)
        for j=1:(S_R1 - S_C1)
            S1 = S1(1:S_C1,1:S_C1);
            C1 = C1(1:S_R1,1:S_C1);
        end
    end 
    if(S_C1>S_R1)
        for j=1:(S_C1 - S_R1)
            S1 = S1(1:S_R1,1:S_R1);
            D1 = D1(1:S_C1,1:S_R1);
        end
    end
    
    
    u = D1*inv(S1)*C1'*Y;     % Residues which gives amplitude and phase
    
    %--------------------


    for j=1:p
            zr(j)=log(Z(j));
            f(j)=[imag(zr(j))/(2*pi*t)];
            d(j)=[real(zr(j))/t];
            sig(j) = d(j);
            w(j) = f(j);
            for l=1:Cit
                AMP(j,l) = 2*abs(u(j,l));   % 2*abs... to count for +ve and -ve freq
                Ph(j,l) = angle(u(j,l));
            end
    end
    
    x = 1;    % counter for adding parameters to excel file
    x2 = 1;   % counter for adding roots and residues
    
   
    for j=r+1:r+p
            Q12(j,1) = x1+1;
            Q12(j,2) = sig(x);
            Q12(j,3) = w(x);
            cc=3;
            for l=1:Cit
                Q12(j,cc+1) = AMP(x,l);
                Q12(j,cc+2) = Ph(x,l)*180/pi;
                Q12(j,cc+3) = AMP(x,l)*exp(sig(x)/w(x));
                cc = cc+3;
            end
            % Factor into consideration 
            % which is D = A*e^(sigma/f) for grouping
                                                
            x=x+1;
            
            M(j,1) = x1+1;
            M(j,2) = z(x2);
            for l=1:Cit
                M(j,2+l) = u(x2,l);
            end
            x2 = x2+1;
    end

    r=r+p;         % updating rows by no of roots

    k = k + N;     % change this to change windowing nature -----------(1)
    
    x1 = x1+1;     % updating window no by 1 only
  end

end

Q11 = {'Window no','Damping','Frequency','Amp1','Phase1','Factor1','Amp2','Phase2','Factor2'};



%xlswrite('out_win1',Q11,'Complete','A1');
%xlswrite('out_win1',Q12,'Complete','A2');   

k2=1;
k=0;
for k3=1:x1
    for j=k+1:k+p
        if (Q12(j,3)>0 && Q12(j,3)<f1)
            for j2=1:Cit
                if (Q12(j,1+3*j2)<0.000001*max(Q12(k+1:k+p,1+3*j2)))
                    Q12(j,1+3*j2) = 0;
                end
            end
            for k1=1:(3*(Cit+1))
                Q13(k2,k1) = Q12(j,k1);
            end
            k2 = k2+1;
        end
    end
    k = k+p;
end

Q13_size = size(Q13,1);

for j2=1:x1
    Q13_count(j2) = 0;
end

for j2=1:x1 
    k=0;
    for k3=1:Q13_size
        if (Q13(k3,1)==j2)
            Q13_count(j2)= Q13_count(j2)+1;
        end
    end
end


k=0;
for k3=1:x1
    for j2=1:Cit
        [Qm(k3,j2),Qn(k3,j2)] = max(Q13((k+1:k+(Q13_count(k3))),3+3*j2));
    end
    [Qmm(k3),Qmn(k3)] = max(Qm(k3));
    k = k+Q13_count(k3);
end




k=0;
for k3=1:x1
    for j2=1+k:Q13_count(k3)+k
        R4(1+k:Q13_count(k3)+k,:) = sortrows(Q13(k+1:k+(Q13_count(k3)),:),-(3+3*Qmn(k3)));
    end
    k = k+Q13_count(k3);
end

k=0;
for k3=1:x1
    for j2=1+k:Q13_count(k3)+k
        for j3=1:3+(4*Cit)
            R5(1+k:Q13_count(k3)+k,j3) = 0;
        end
        for j3=1:3+(3*Cit)
            R5(1+k:Q13_count(k3)+k,j3) = R4(1+k:Q13_count(k3)+k,j3);
        end
        for j3=1:Cit
            R5(1+k:Q13_count(k3)+k, j3+3+(3*Cit)) = R4(1+k:Q13_count(k3)+k,2+3*j3) - R4(1+k:Q13_count(k3)+k,2+3*Qmn(k3));
        end
    end
    k = k+Q13_count(k3);
end
        
        
R6 = {'Window no','Damping','Frequency','Amp1','Phase1','Factor1','Amp2','Phase2','Factor2','PhDiff1','PhDiff2'};        

xlswrite('sorted_multi',R6,'Complete','A1');
xlswrite('sorted_multi',R5,'Complete','A2');

for j=1:W/N
    k=0;
    for j1 = 1:size(R4,1)
        if(R4(j1,1)==j)
            k = k+1;
        end
    end
    R8(j)=k;
end

rc = 1;
for j=1:W/N
R9(j,1) = R4(rc,3);
%R9(j,1) = R4(rc,3);
R9(j,2) = R4(rc,6);
% R9(j,3) = R4(rc,9);
% R9(j,4) = R4(rc,12);
rc = rc + R8(j);
end

for j=1:W/N
    R10(j) = j;
end

scatter(R10,R9(:,1));
grid on
