clear all
clc
 %Semester project for solving APc1-5; A 2-dimensional Poisson equation
 tic
 
 
 
 Nodes=200;
 x = length(Nodes);
 y= length(Nodes);
 dx= (pi+pi)/(Nodes+1);
 dy=dx;
 %defining step intervals 
 x= -pi:dx:pi;
 y=-pi:dy:pi;
 %no need to compute this every time in loop when doing gauss-siedel or
 %over relaxation
 dx2=dx*dx;
%  k=-.01
%  dx2=dx2/k;
 %use vectorization here
 boundary_PHI= 1*((y-(-pi)).*(y-(-pi))).*sin( (pi/2) * ((y-(-pi))/(pi-(-pi)))  );  % goes on U(1,:)
 boundary_PSI= 1*cos(pi.*(y-(-pi))).*cosh(pi-y);%goes on U(Nodes+2,:)
 
 F= zeros(Nodes+2);
 for j=2:Nodes+1
     for i=2:Nodes+1
         F(i,j)= -1*sin((x(i)+pi)/(2) ) *cos(pi*.5*( ((y(j)+pi)/(pi)) +1));
     end
 end
 %put zero at 'edges'
%  F(1,:)=0;
%  F(Nodes+2,:)=0;
%  F(:,1)=0;
%  F(:,Nodes+2)=0;
 U_old = zeros(Nodes+2);
 U = ones(Nodes+2);
 %U = zeros(Nodes+2);
 % put on BCs for 'padding' 
 U(1,:)=boundary_PHI;
 U(Nodes+2,:)=boundary_PSI;
 Comparer=zeros(Nodes+2,Nodes);
 counter_gauss=0;
 condition_gauss=0;
 
% ******** beginning of Gauss-Seidel ********

%  while ( condition_gauss ==0) %Iteration number 
%        for i=2:Nodes+1%must solve for U in here as Ditchelt B.C not specified- U will be unknown here
%          % moved to solve for boundary u here- had if statements and took too long
%            U(i,1)= ( U(i-1,1) + 2*U(i,2) + U(i+1,1) - dx2*F(i,1) )/4; 
%         %Managed to cut test time to 663 sec from 1500 sec for very large
%         %node quantity
%        end
%        
%        for j=2:Nodes+1
%         for i=2:Nodes+1 % U at boundary given, dont have to solve for them
%            U(i,j)= ( U(i-1,j) + U(i,j+1) + U(i+1,j) + U(i,j-1) - dx2*F(i,j) )/4;
%                
%         end
%        end
%        
%        for i=2:Nodes+1 %must solve for U in here as Ditchelt B.C not specified- U will be unknown here
%              
%            U(i,Nodes+2)= ( U(i-1,Nodes+2) + U(i+1,Nodes+2) + 2*U(i,Nodes+1)- dx2*F(i,j) )/4; %made two additional for loops rather than having
%              % only one main for loop for i & a bunch of if statements that
%              % had to be checked every time
%        
%        end
%        for j=1:Nodes+2
%            for i=2:Nodes+1 
%                
%            Comparer(i,j)= abs( ( U(i,j)-U_old(i,j) ) / U(i,j));
%            
%            end
%        end
%        if max(Comparer) <.001
%            
%            condition_gauss=1; % Changes condition to 1 to get out of while loop
%       
%        else
%            
%            U_old=U;
%            counter_gauss=counter_gauss+1;
%        
%        end
%                
%            
%  end
%  counter_gauss
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % ******** end of Gauss-Seidel solver ******** %
 




%******** Beginning of SOR solver ********

% condition_SOR =0;
% counter_SOR=0;
% beta= 1.5; % recommended value for successive over relaxation
% alpha= beta/4; % put at multiplier to save some time 
% gamma = 1-beta;
% %   for k=1:100;
%  while ( condition_SOR ==0 ) %Iteration number 
%        for i=2:Nodes+1%must solve for U in here as Ditchelt B.C not specified- U will be unknown here
%          
%            U(i,1)= U(i,1)*gamma + ( U(i-1,1) + 2*U(i,2) + U(i+1,1) -dx2*F(i,1) )* alpha; 
%        
%        end
%        
%        for j=2:Nodes+1
%         for i=2:Nodes+1 % U at boundary given, dont have to solve for them
%            U(i,j)= U(i,j)*gamma +  ( U(i-1,j) + U(i,j+1) + U(i+1,j) + U(i,j-1) -dx2*F(i,j) )* alpha;
%                
%         end
%        end
%        
%        for i=2:Nodes+1 %must solve for U in here as Ditchelt B.C not specified- U will be unknown here
%              
%            U(i,Nodes+2)= U(i,Nodes+2)*gamma+ (  U(i-1,Nodes+2) + U(i+1,Nodes+2) + 2*U(i,Nodes+1) - dx2*F(i,j) )*alpha; %made two additional for loops rather than having
%              % only one main for loop for i & a bunch of if statements that
%              % had to be checked every time
%        
%        end
%        for j=1:Nodes+2
%            for i=2:Nodes+1 
%                
%            Comparer(i,j)= abs( ( U(i,j)-U_old(i,j) )) / abs(U(i,j));
%            
%            end
%        end
%        if max(Comparer) <.0001
%            
%            condition_SOR=1; % Changes condition to 1 to get out of while loop
%       
%        else
%            
%            U_old=U;
%            counter_SOR=counter_SOR+1;
%        
%        end
%                
%            
%  end
% counter_SOR
 %******** End of SOR SOLVER ********
 
% U_Plot = transpose(U);
% surf(U_Plot);
%  title('Surface Plot for solution of given Poisson equation')
%  ylabel('y-axis')
%  xlabel('x-axis')
%  

 
% plot(x,U(1,:))

%Method of Manufactured Solutions Comparison
% From Professor's notes, use the following:
% V will be subsitituted in to 'manufacture' a solution
% Q will be laplacian of manufactured solution
V=zeros(Nodes+2,Nodes+2);
U= zeros(Nodes+2,Nodes+2);
Q= zeros(Nodes+2,Nodes+2);
for j=1: Nodes+2
    for i=1:Nodes+2
        V(i,j)= ((x(i))^3)*((2*pi-x(i))^3) * ((y(j))^3)*((2*pi-y(j))^3); % solving for made up V on all points on 'plate'
        % Q will be laplacianV which will be used to help solve for answer
        % numerically
          
    end 
end
for j=1:Nodes+2
    for i=2:Nodes+1
       Q(i,j) = ( 6*x(i)*( 2*pi-x(i))^3 + 6*(x(i)^3)*(2*pi - x(i)) -18*(x(i)^2)*(2*pi-x(i))^2)*((y(j))^3)*(2*pi - y(j))^3  + ( 6*y(j)*(2*pi-y(j))^3 + 6*((y(j))^3)*(2*pi-y(j))-18*(y(j)^2)*(2*pi-y(j))^2 ) * (x(i)^3)*(2*pi-x(i))^3; 
    end
end
% at Dirichlet boundaries, values given 
U(1,:)= V(1,:);
U(Nodes+2,:)= V(Nodes+2,:);
    Moon=((3*(-pi)^2)*((2*pi+pi)^3)   -3*((-pi)^3)*((2*pi+pi)^2)); %solve for at y=-pi
    Sun = (  (3*(pi^2)*((pi)^3)   -3*(pi^3)*((pi)^2))  ); %solves for at y=pi;
        dV_zero=     Moon.*((2*pi-x).^3) .* (x .^3); 
        dV_N_PLUS=   Sun .*((2*pi-x).^3 .*(x .^3)); 

% since now dv/dy not equal to zero, will need value of dx*2....

two_dx= 2*dx;

% ******** beginning of Gauss-Seidel To test MMS ********

 while ( condition_gauss ==0) %Iteration number 
       for i=2:Nodes+1%must solve for U in here as Ditchelt B.C not specified- U will be unknown here
         
          U(i,1)= ( U(i-1,1) + 2*U(i,2)-dV_zero(i)*two_dx + U(i+1,1) - dx2*Q(i,1) )/4; 
        
       end
       
       for j=2:Nodes+1
        for i=2:Nodes+1 % U at boundary given, dont have to solve for them
           U(i,j)= ( U(i-1,j) + U(i,j+1) + U(i+1,j) + U(i,j-1) - dx2*Q(i,j) )/4;
               
        end
       end
       
       for i=2:Nodes+1 %must solve for U in here as Ditchelt B.C not specified- U will be unknown here
             
           U(i,Nodes+2)= ( U(i-1,Nodes+2) + U(i+1,Nodes+2) +dV_N_PLUS(i)*two_dx + 2*U(i,Nodes+1) - dx2*Q(i,j) )/4; %made two additional for loops rather than having
             % only one main for loop for i & a bunch of if statements that
             % had to be checked every time
       
       end
       for j=1:Nodes+2
           for i=1:Nodes+2 
               
           Comparer(i,j)= abs( ( U(i,j)-U_old(i,j) ) /  U(i,j));
           
           end
       end
       if max(Comparer) <.0000000001
           
           condition_gauss=1; % Changes condition to 1 to get out of while loop
      
       else
           
           U_old=U;
           counter_gauss=counter_gauss+1;
       
       end
               
           
 end
 counter_gauss
U_Plot = transpose(U);
V_Plot = transpose(V);
for i=1:Nodes+2
    for j=1:Nodes+2
        comparer_2(i,j)= abs((U_Plot(i,j)-V_Plot(i,j))/V_Plot(i,j));
    end
end
Error = V_Plot - U_Plot;
 surf(U_Plot) 
hold on
 surf(V_Plot) 
 title('Surface Plot for solution of given Poisson equation')
 ylabel('y-axis')
 xlabel('x-axis')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      toc
 