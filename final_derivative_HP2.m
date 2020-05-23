

function derivative_xi=final_derivative_HP2(t,y,n,f,initial_delta_x,av_nb_of_nodes,k,typical_surf)
   

   
    #---values of Xi derivative/stress for first and last cells
   derivative_xi(1,1)=(k*((y(2)-y(1)-initial_delta_x)/initial_delta_x)-f)/n;
   derivative_xi(av_nb_of_nodes,1)=(-k*((y(av_nb_of_nodes)-y(av_nb_of_nodes-1)-initial_delta_x)/initial_delta_x)+f)/n;
   #---values of Xi derivative for other cells, self propelling force depending on position
    
   if mod(av_nb_of_nodes,2)==0
     
     for i = 2:round(av_nb_of_nodes/2)
        derivative_xi(i,1)=(-k/initial_delta_x*(y(i)-y(i-1)-initial_delta_x)+k/initial_delta_x*(y(i+1)-y(i)-initial_delta_x)-f)/n;
     endfor
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        derivative_xi(i,1)=(-k/initial_delta_x*(y(i)-y(i-1)-initial_delta_x)+k/initial_delta_x*(y(i+1)-y(i)-initial_delta_x)+f)/n;
     endfor
     derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-k/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
   
   endif
   
   
   if mod(av_nb_of_nodes,2)==1
    
    for i = 2:round(av_nb_of_nodes/2)-1
        derivative_xi(i,1)=(-k/initial_delta_x*(y(i)-y(i-1)-initial_delta_x)+k/initial_delta_x*(y(i+1)-y(i)-initial_delta_x)-f)/n;
     endfor
     for i = round(av_nb_of_nodes/2)+2:av_nb_of_nodes-1
        derivative_xi(i,1)=(-k/initial_delta_x*(y(i)-y(i-1)-initial_delta_x)+k/initial_delta_x*(y(i+1)-y(i)-initial_delta_x)+f)/n;
     endfor
     derivative_xi(round(av_nb_of_nodes/2)+1,1)=(-k/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x)+k/initial_delta_x*(y(round(av_nb_of_nodes/2)+2)-y(round(av_nb_of_nodes/2)+1)-initial_delta_x))/n;  
     derivative_xi(round(av_nb_of_nodes/2),1)=(-k/initial_delta_x*(y(round(av_nb_of_nodes/2))-y(round(av_nb_of_nodes/2-1))-initial_delta_x)+k/initial_delta_x*(y(round(av_nb_of_nodes/2)+1)-y(round(av_nb_of_nodes/2))-initial_delta_x))/n;  
  
   endif
   
   
   derivative_xi(av_nb_of_nodes,1)=(-k/initial_delta_x*(y(av_nb_of_nodes)-y(av_nb_of_nodes-1)-initial_delta_x)+f)/n;


   
end