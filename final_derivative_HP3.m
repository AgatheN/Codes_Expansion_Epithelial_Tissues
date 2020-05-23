
function derivative_xi=final_derivative_HP3(t,y,n,f,initial_delta_x,av_nb_of_nodes,k,typical_surf,epsilon_s,flag_threshold_reached,T_eps_s)
   
   
   global flag_threshold_reached;
   global T_eps_s;
   
   derivative_xi(1,1)=((k/initial_delta_x)*(y(2)-y(1)-initial_delta_x)-f)/n;
   derivative_xi(av_nb_of_nodes,1)=((-k/initial_delta_x)*(y(av_nb_of_nodes)-y(av_nb_of_nodes-1)-initial_delta_x)+f)/n;
   
   
   #---values of Xi derivative, self propelling force depending on strain of neighbouring cells
   if mod(av_nb_of_nodes,2)==0  #pair nb of nodes
     for i = 2:round(av_nb_of_nodes/2) #before the midline 
         
         if flag_threshold_reached(i)==1   #if strain threshold is already reached
           derivative_xi(i,1)=(-f+(-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
          #self propelling force is exerced
         else
              if ((y(i)-y(i-1)-initial_delta_x)/(initial_delta_x))>=epsilon_s
                  T_eps_s(i)=t;
                  flag_threshold_reached(i)=1;#if strain threshold is already reached
                  derivative_xi(i,1)=(-f+(-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
                  #self propelling force is exerced
              else
                  derivative_xi(i,1)=((-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
                  #self propelling force is not exerced
             endif
         endif
     endfor
     
     
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1 #after the midline
       #same behaviour as described before
        if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(f-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n; 
        else
      
            if(y(i+1)-y(i)-initial_delta_x)/(initial_delta_x)>=epsilon_s 
                T_eps_s(i)=t;
                flag_threshold_reached(i)=1;
                derivative_xi(i,1)=(f-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            else
                derivative_xi(i,1)=(-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            endif
        endif
     endfor
     
      derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-(k/initial_delta_x)*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+(k/initial_delta_x)*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
   endif
   
   #if odd nb of nodes, same behaviour
   if mod(av_nb_of_nodes,2)==1
     for i = 2:round(av_nb_of_nodes/2)-1
         if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(-f+(-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
         else
              if ((y(i)-y(i-1)-initial_delta_x)/(initial_delta_x))>=epsilon_s
                 T_eps_s(i)=t;
               
                  flag_threshold_reached(i)=1;
                  derivative_xi(i,1)=(-f+(-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              else
                  derivative_xi(i,1)=((-k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
              endif
         endif
     endfor
     
     
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        if flag_threshold_reached(i)==1
           derivative_xi(i,1)=(f-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n; 
        else
            if(y(i+1)-y(i)-initial_delta_x)/(initial_delta_x)>=epsilon_s 
               T_eps_s(i)=t;
               
                flag_threshold_reached(i)=1;
                derivative_xi(i,1)=(f-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            else
                derivative_xi(i,1)=(-(k/initial_delta_x)*(y(i)-y(i-1)-initial_delta_x)+(k/initial_delta_x)*(y(i+1)-y(i)-initial_delta_x))/n;
            endif
        endif
     endfor
     derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-k/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
     derivative_xi(round(av_nb_of_nodes/2),1)=(-k/initial_delta_x*(y(round(av_nb_of_nodes/2))-y(round(av_nb_of_nodes/2-1))-initial_delta_x)+k/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x))/n;  
   endif

end
