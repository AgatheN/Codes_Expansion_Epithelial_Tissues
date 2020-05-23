function derivative_xi=final_derivative_HP1(t,y,n,f,initial_delta_x,av_nb_of_nodes,k)


  #---values of Xi derivative/stress for first and last cells
  derivative_xi(1,1)=(k*((y(2)-y(1)-initial_delta_x)/initial_delta_x)-f)/n;
  derivative_xi(av_nb_of_nodes,1)=(-k*((y(av_nb_of_nodes)-y(av_nb_of_nodes-1)-initial_delta_x)/initial_delta_x)+f)/n;
    
  #---values of Xi derivative/strain rate/stress for other cells
  for i = 2:av_nb_of_nodes-1
      derivative_xi(i,1)=(-k*((y(i)-y(i-1)-initial_delta_x)/initial_delta_x)+k*((y(i+1)-y(i)-initial_delta_x)/initial_delta_x))/n;
      endfor 
  
end