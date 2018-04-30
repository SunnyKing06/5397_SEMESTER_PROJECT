clear all
clc
 %Semester project for solving APc1-5; A 2-dimensional Poisson equation
 tic
 Nodes= 400;
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
 
 %use vectorization here
 boundary_PHI= 1*((y-(-pi)).*(y-(-pi))).*sin( (pi/2) * ((y-(-pi))/(pi-(-pi)))  );  % goes on U(1,:)
 boundary_PSI= 1*cos(pi.*(y-(-pi))).*cosh(pi-y);%goes on U(Nodes+2,:)
 
 F= zeros(Nodes+2);
 for j=1:Nodes+2
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
 counter=0;
 condition=0;
 
 % beginning of Gauss-Seidel
%  for k=1:200;
 while ( condition ==0) %Iteration number 
       for i=2:Nodes+1%must solve for U in here as Ditchelt B.C not specified- U will be unknown here
         % moved to solve for boundary u here- had if statements and took too long
           U(i,1)= ( U(i-1,1) + 2*U(i,2) + U(i+1,1) + -dx2*F(i,1) )/4; 
        %Managed to cut test time to 663 sec from 1500 sec for very large
        %node quantity
       end
       
       for j=2:Nodes+1
        for i=2:Nodes+1 % U at boundary given, dont have to solve for them
           U(i,j)= ( U(i-1,j) + U(i,j+1) + U(i+1,j) + U(i,j-1) -dx2*F(i,j) )/4;
               
        end
       end
       
       for i=2:Nodes+1 %must solve for U in here as Ditchelt B.C not specified- U will be unknown here
             
           U(i,Nodes+2)= ( U(i-1,Nodes+2) + U(i+1,Nodes+2) + 2*U(i,Nodes+1)+ -dx2*F(i,j) )/4; %made two additional for loops rather than having
             % only one main for loop for i & a bunch of if statements that
             % had to be checked every time
       
       end
       for j=2:Nodes+1
           for i=1:Nodes+2 
               
           Comparer(i,j)= abs( ( U(i,j)-U_old(i,j) )) / abs(U(i,j));
           
           end
       end
       if max(Comparer) <.001
           
           condition=1; % Changes condition to 1 to get out of while loop
      
       else
           
           U_old=U;
           counter=counter+1;
       
       end
               
           
 end
 counter
 
 
% end of Gauss-Seidel solver
 


%Beginning of Jacobi linear solver

% U_old=U; % Old starting as assumed guess of U being zero plus boundary condition
% for k=1:100  %Iteration number 
%     for j=1:Nodes+2 %must solve for U in here as Ditchelt B.C not specified- U will be unknown here
%         for i=2:Nodes+1 % U at boundary given, dont have to solve for them
%             if (j==1) || (j==(Nodes+2)) % check to see if at boundary
%                 if j==1
%                 U(i,j)= ( U_old(i-1,j) + 2*U_old(i,j+1) + U_old(i+1,j) +  -dx2*F(i,j) )/4;
%                 else
%                 U(i,j)= ( U_old(i-1,j) + U_old(i+1,j) + 2*U_old(i,j-1) -dx2*F(i,j) )/4;
%                 end
%             else
%           U(i,j)= ( U_old(i-1,j) + U_old(i,j+1) + U_old(i+1,j) + U_old(i,j-1) -dx2*F(i,j) )/4;
%             end
%         end
%     end
%     U_old=U; % new becomes old at end of iteration, can use on the next turn to get new ones
%  end

% % end of Jacobi linear solver

 U_Plot = transpose(U);
 
 surf(U_Plot)
 ylabel('y-axis')
 xlabel('x-axis')
toc
 
% plot(x,U(1,:))
     
 