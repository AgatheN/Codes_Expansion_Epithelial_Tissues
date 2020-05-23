
function derivative_xi=final_derivative_HP4(t,y,n,f,initial_delta_x,av_nb_of_nodes,av_nb_of_cells,k,epsilon_s,K,delta_t,epsilon_th,t_threshold,polynome,t_threshold_matrix,flag_threshold_reached,T_eps_s,count_t)



   global K;
   global t_threshold;
   global t_threshold_matrix;
   global flag_threshold_reached;
   global T_eps_s;
   global count_t;


  

for i=1:av_nb_of_cells
    if 0<=(t-t_threshold(i)) && (t-t_threshold(i))<=delta_t 
       #If the reinforcement/fluidization cycle is not finished yet
        K(i)= polyval(polynome,(t-t_threshold(i))); #value of the cycle at time t
    else  
         if ((y(i+1)-y(i)-initial_delta_x)/initial_delta_x)>epsilon_th 
         #If strain threshold eps_th is reached for the given cell i
               t_threshold(i)=t; #t is defined as the time at which the threshold is reached
               K(i)= polyval(polynome,(t-t_threshold(i)));
               t_threshold_matrix(count_t,i)=t; 
               #t is stored to further get the values of V,S,E
         else 
              K(i)=k; #k is equal to its standard value
         endif
    endif
endfor
  
  

derivative_xi(1,1)=((K(1)/initial_delta_x)*(y(2)-y(1)-initial_delta_x)-f)/n;
derivative_xi(av_nb_of_nodes,1)=((-K(av_nb_of_cells)/initial_delta_x)*(y(av_nb_of_nodes)-y(av_nb_of_nodes-1)-initial_delta_x)+f)/n;
   
   
   #---values of Xi derivative, self propelling force depending on strain of neighbouring cells
   if mod(av_nb_of_nodes,2)==0
     for i = 2:round(av_nb_of_nodes/2)
         
         if flag_threshold_reached(i)==1 
         
           derivative_xi(i,1)=(-f+(-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
         else
              if ((y(i)-y(i-1)-initial_delta_x)/(initial_delta_x))>=epsilon_s
                  T_eps_s(i)=t;
                  flag_threshold_reached(i)=1;
                  derivative_xi(i,1)=(-f+(-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              else
                  derivative_xi(i,1)=((-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              endif
         endif
     endfor
     
     
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(f-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n; 
        else
      
            if(y(i+1)-y(i)-initial_delta_x)/(initial_delta_x)>=epsilon_s 
                T_eps_s(i)=t;
                flag_threshold_reached(i)=1;
                derivative_xi(i,1)=(f-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            else
                derivative_xi(i,1)=(-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            endif
        endif
     endfor
     
      derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-(K(round(av_nb_of_nodes/2))/initial_delta_x)*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+(K(round(av_nb_of_nodes/2)+1)/initial_delta_x)*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
   endif
   
   if mod(av_nb_of_nodes,2)==1
     for i = 2:round(av_nb_of_nodes/2)-1
         if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(-f+(-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
         else
              if ((y(i)-y(i-1)-initial_delta_x)/(initial_delta_x))>=epsilon_s
                 T_eps_s(i)=t;
               
                  flag_threshold_reached(i)=1;
                  derivative_xi(i,1)=(-f+(-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              else
                  derivative_xi(i,1)=((-K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              endif
         endif
     endfor
     
     
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(f-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n; 
        else
            if(y(i+1)-y(i)-initial_delta_x)/(initial_delta_x)>=epsilon_s 
               T_eps_s(i)=t;
               
                flag_threshold_reached(i)=1;
                derivative_xi(i,1)=(f-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            else
                derivative_xi(i,1)=(-(K(i-1)/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(K(i)/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            endif
        endif
     endfor
     derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-K(round(av_nb_of_nodes/2))/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+K(round(av_nb_of_nodes/2)+1)/initial_delta_x*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
     derivative_xi(round(av_nb_of_nodes/2),1)=(-K(round(av_nb_of_nodes/2)-1)/initial_delta_x*(y(round(av_nb_of_nodes/2))-y(round(av_nb_of_nodes/2-1))-initial_delta_x)+K(round(av_nb_of_nodes/2))/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x))/n;  
   endif
   
count_t+=1;
end