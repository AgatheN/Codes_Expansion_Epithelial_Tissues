

function [V,S,E_point]=Calc_v_s_e_HP4(T,Xi,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,typical_surf,epsilon_s,T_eps_s,t_threshold_matrix,polynome,delta_t,epsilon_th,t_max)

#Initialisations of matrices / arrays

V=zeros(size(T),av_nb_of_nodes);
S=zeros(size(T),av_nb_of_cells); 
E_point=zeros(size(T),av_nb_of_cells);
t_threshold_eps_t=ones(av_nb_of_cells,1)*t_max
flag_threshold_reached=zeros(av_nb_of_nodes,1);
   

for t =1:size(T) 
    #Values of k for each cell, time T(t)
    for i=1:av_nb_of_cells
      if t==1 #initialisation of K array
         K=ones(av_nb_of_cells,1)*k;
      else  
          if 0<=(T(t)-t_threshold_eps_t(i)) && (T(t)-t_threshold_eps_t(i))<=delta_t 
               #If the reinforcement/fluidization cycle is not finished yet
               K(i)= polyval(polynome,(T(t)-t_threshold_eps_t(i))); #value of the cycle at time T(t)
          else  
               indx_t_threshold_i=find(T(t-1)<t_threshold_matrix(:,i) & t_threshold_matrix(:,i)<=T(t),1); #&& doesn't wor
               # time index of the matrix giving the times at which a cycle starts for each cell
               # index isn't empty if there the cell starts a cycle between T(t-1) and T(t)
               if isempty(indx_t_threshold_i)  #If the cell doesn't start a reinforcement cycle
                   K(i)=k; #k is equal to its standard value
               else
                   t_threshold_eps_t(i)=t_threshold_matrix(indx_t_threshold_i,i);
                   #Exact time between T(t-1) and T(t) when the cycle begins
                   K(i)= polyval(polynome,(T(t)-t_threshold_eps_t(i)));   
                   #Evaluation of k thanks to the polynomic function
               endif
          endif
       endif
    endfor
     
    
     
     
     V(t,1)=((K(1)/initial_delta_x)*(Xi(t,2)-Xi(t,1)-initial_delta_x)-f)/n;
     V(t,av_nb_of_nodes)=((-K(av_nb_of_cells)/initial_delta_x)*(Xi(t,av_nb_of_nodes)-Xi(t,av_nb_of_nodes-1)-initial_delta_x)+f)/n;
     S(t,1)=K(1)*((Xi(t,2)-Xi(t,1)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
     
     
     #---values of Xi derivative, self propelling force depending on strain of neighbouring cells
  
  #FOR EXPLANATIONS SEE Calc_v_s_e_HP3 : self propelling forces is exerced if a strain threshold epsilon_s is reached
  
  # S is modified through this iteration loop since it depends on k
  
  if mod(av_nb_of_nodes,2)==0
    
     for i = 2:round(av_nb_of_nodes/2) 
          if flag_threshold_reached(i)==1
             V(t,i)=(-f+(-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
          elseif flag_threshold_reached(i)==0
              if T(t)>=T_eps_s(i) 
                  flag_threshold_reached(i)=1;
                  V(t,i)=(-f+(-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              else
                  V(t,i)=((-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              endif
         endif
         S(t,i)=K(i)*((Xi(t,i+1)-Xi(t,i)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
     endfor
       
       
     for i = round(av_nb_of_nodes/2)+1:av_nb_of_nodes-1
       if flag_threshold_reached(i)==1
           V(t,i)=(f-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n; 
       elseif flag_threshold_reached(i)==0
            if T(t)>=T_eps_s(i) 
                flag_threshold_reached(i)=1;
                V(t,i)=(f-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            else
                V(t,i)=(-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            endif
        endif
        S(t,i)=K(i)*((Xi(t,i+1)-Xi(t,i)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
     endfor
       
       V(t,round(av_nb_of_nodes/2)+1)=(-(K(round(av_nb_of_nodes/2))/initial_delta_x)*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+(K(round(av_nb_of_nodes/2)+1)/initial_delta_x)*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
       
 endif


 
 if mod(av_nb_of_nodes,2)==1
       for i = 2:round(av_nb_of_nodes/2)-1
          if flag_threshold_reached(i)==1
           V(t,i)=(-f +(-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
          elseif flag_threshold_reached(i)==0
              if T(t)>=T_eps_s(i) 
                  flag_threshold_reached(i)=1;
                  V(t,i)=(-f+(-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              else
                  V(t,i)=((-K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
              endif
         endif
         S(t,i)=K(i)*((Xi(t,i+1)-Xi(t,i)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
     endfor
       
       
     for i = round(av_nb_of_nodes/2)+1:av_nb_of_nodes-1
       if flag_threshold_reached(i)==1
           V(t,i)=(f-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n; 
        elseif flag_threshold_reached(i)==0
            if T(t)>=T_eps_s(i) 
                flag_threshold_reached(i)=1;
                V(t,i)=(f-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            else
                V(t,i)=(-(K(i-1)/initial_delta_x)*(Xi(t,i)-Xi(t,i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(Xi(t,i+1)-Xi(t,i)-initial_delta_x))/n;
            endif
        endif
        S(t,i)=K(i)*((Xi(t,i+1)-Xi(t,i)-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);
     endfor
     V(t,round(av_nb_of_nodes/2)+1)=(-K(round(av_nb_of_nodes/2))/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)+K(round(av_nb_of_nodes/2)+1)/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+2)-Xi(t,round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
     V(t,round(av_nb_of_nodes/2))=(-K(round(av_nb_of_nodes/2)-1)/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2))-Xi(t,round(av_nb_of_nodes/2-1))-initial_delta_x)+K(round(av_nb_of_nodes/2))/initial_delta_x*(Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x))/n;  
     S(t,round(av_nb_of_nodes/2))=K(round(av_nb_of_nodes/2))*((Xi(t,round(av_nb_of_nodes/2)+1)-Xi(t,round(av_nb_of_nodes/2))-initial_delta_x)/initial_delta_x)*10^-9/(typical_surf*10^-12);

  endif
     
    
     
 #-----values of strain rate/stress for other cells
 for i = 2:av_nb_of_nodes-1
       E_point(t,i-1)=(V(t,i)-V(t,i-1))/(initial_delta_x*60);
 endfor
    
 #-----value of strain rate for the last cell
 E_point(t,av_nb_of_nodes-1)=(V(t,av_nb_of_nodes)-V(t,av_nb_of_nodes-1))/(initial_delta_x*60);
  
  endfor
end
