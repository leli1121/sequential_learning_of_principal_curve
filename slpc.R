library(numDeriv)

slpc = function(vertices, past_obs, new_obs, p, R, L, delta, c0, c1, c2, c3, n_nearest = 3,  learning_rate = 0.5, epsilon = 0.5, lat){
  
  whole_lattice = lattice_generator_R(R, delta, d)
  
  
  if (class(vertices) != 'matrix'){
    print ('vertices should have matrix form !')
  }
  
  nbr_vertices  = dim(vertices)[1]
  
  # When nbr_vertices >=3, k_neighbors takes value in [nbr_vertices-1, nbr_vertices, nbr_vertices+1] and it would decide whether to delete, keep
  # or add one vertices.
  k_neighbors = seq(max(2, (nbr_vertices-1)), min((nbr_vertices +1), p))
  R_local = 2*(sqrt(max(rowSums(rbind(past_obs, new_obs)^2))) + 2)/sqrt(d)
  # round t
  t_length = dim(rbind(past_obs, new_obs))[1]
  inverse_eta_t = c0*d*(2*R_local+delta)^2*sqrt((exp(1)-1)*t_length)/sqrt(h(vertices, L, c1, c2, c3, p, R_local, delta, indicator = TRUE)) 
  k_final_vertices = list()
  
  for (k in k_neighbors){
    new_k_vertices = add_keep_delete_vertex(vertices, k, new_obs)
    new_k_nbr = dim(new_k_vertices)[1]
    
    V_1_names = sapply(1:new_k_nbr, function(x) paste('vertex_', x, sep='')) #unlist(lapply(1:new_k_nbr, function(x) paste('vertex_', x, sep='')))
    V_1_segments_names = sapply(1:(new_k_nbr-1), function(x) paste('segment_', x, sep=''))
    
    if (new_k_nbr >=5){
      case = 'vtc_5'
    }else if ((2 <=new_k_nbr) & (new_k_nbr<=4)){
      case = paste('vtc_', new_k_nbr, sep = '')
    }else{
      case = 'do not exist !'
    }
    
    switch(case,
           vtc_5 = {
             partition_result = voronoi_partitions(new_k_vertices, past_obs, new_obs)
             voronoi_of_all_obs = partition_result$partition_of_all_observations 
             index_of_new_obs = partition_result$index_of_obs
             neighbor_names = partition_result$neighbor_names
             local_obs_list =  lapply(1:length(neighbor_names), function(i) voronoi_of_all_obs[neighbor_names[i]][[1]])
             
             
             if (index_of_new_obs %in% c(V_1_names[1], V_1_segments_names[1])){# index_of_new_obs is either vertex_1 or segement_1
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex or segment
               ith = strtoi(name_and_index[2]) ## ith = 1
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               update_v1_v2 = function(local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
             }else if(index_of_new_obs == V_1_names[2]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 2
               initial_local_v = new_k_vertices[1:3, ]
               
               update_v1_v2_v3 = function(local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
               }
               
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
             }else if(index_of_new_obs %in% V_1_segments_names[2:(new_k_nbr-2)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith:ith+1
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 2 -ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-1), ], end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
               
             }else if(index_of_new_obs %in% V_1_names[3:(new_k_nbr-2)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith-1 : ith+1
               initial_local_v = new_k_vertices[(ith-1):(ith+1), ]
               
               update_3_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
             }else if(index_of_new_obs == V_1_names[(new_k_nbr-1)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## k-1 : k+1
               
               initial_local_v = new_k_vertices[(ith-1):(ith+1), ]
               
               update_3_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               
               Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)           
               
             }else if(index_of_new_obs == V_1_segments_names[new_k_nbr-1]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith = k, k : k+1
               
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-1), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)             
               
             }else{
               ith = new_k_nbr
               initial_local_v = new_k_vertices[(ith-1):ith, ]
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)           
               
             }
           },
           
           
           vtc_4 = {
             partition_result = voronoi_partitions(new_k_vertices, past_obs, new_obs)
             voronoi_of_all_obs = partition_result$partition_of_all_observations 
             index_of_new_obs = partition_result$index_of_obs
             neighbor_names = partition_result$neighbor_names
             local_obs_list =  lapply(1:length(neighbor_names), function(i) voronoi_of_all_obs[neighbor_names[i]][[1]])
             if (index_of_new_obs %in% c(V_1_names[1], V_1_segments_names[1])){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 1
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               #start_v =NULL
               #end_v = new_k_vertices[(ith+2):new_k_nbr,]
               #local_v = new_k_vertices[ith:(ith+1), ]
               update_v1_v2 = function(local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t) 
               
             }else if(index_of_new_obs == V_1_names[2]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 2
               initial_local_v = new_k_vertices[1:3, ]
               
               update_v1_v2_v3 = function(local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
               
             }else if(index_of_new_obs %in% V_1_segments_names[2:(new_k_nbr-2)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith:ith+1
               
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 2 - ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-1), ], end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)           
               
               
               
               
             }else if(index_of_new_obs == V_1_names[(new_k_nbr-1)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## k-1 : k+1
               
               initial_local_v = new_k_vertices[(ith-1):(ith+1), ]
               
               update_3_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               
               Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)           
               
               
             }else if(index_of_new_obs == V_1_segments_names[new_k_nbr-1]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith = k, k : k+1
               
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-1), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)   
               
             }else{
               ith = new_k_nbr
               initial_local_v = new_k_vertices[(ith-1):ith, ]
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)           
               
             }
           },
           
           
           
           vtc_3 = {
             partition_result = voronoi_partitions(new_k_vertices, past_obs, new_obs)
             voronoi_of_all_obs = partition_result$partition_of_all_observations 
             index_of_new_obs = partition_result$index_of_obs
             neighbor_names = partition_result$neighbor_names
             local_obs_list =  lapply(1:length(neighbor_names), function(i) voronoi_of_all_obs[neighbor_names[i]][[1]])
             if (index_of_new_obs %in% c(V_1_names[1], V_1_segments_names[1])){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 1
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               update_v1_v2 = function(local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,]){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = new_k_vertices[(ith+2):new_k_nbr,], whole_lattice, n_nearest, inverse_eta_t)
               
               
             }else if(index_of_new_obs == V_1_names[2]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 2
               initial_local_v = new_k_vertices[1:3, ]
               
               update_v1_v2_v3 = function(local_v, start_v = NULL, end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2_v3, initial_local_v, start_v = NULL, end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
               
               
               
               
             }else if(index_of_new_obs == V_1_names[(new_k_nbr-1)]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## k-1 : k+1
               
               initial_local_v = new_k_vertices[(ith-1):(ith+1), ]
               
               update_3_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               
               Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_3_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
               
             }else if(index_of_new_obs == V_1_segments_names[new_k_nbr-1]){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## 
               ith = strtoi(name_and_index[2]) ## ith = k, k : k+1
               
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 2-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-1), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-1), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
               
               
             }else{
               ith = new_k_nbr
               initial_local_v = new_k_vertices[(ith-1):ith, ]
               update_2_vertices = function(local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2]) + 3-ith
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = new_k_vertices[1:(ith-2), ] , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = new_k_vertices[1:(ith-2), ], end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
             }
           },
           
           
           
           vtc_2 = {
             partition_result = voronoi_partitions(new_k_vertices, past_obs, new_obs)
             voronoi_of_all_obs = partition_result$partition_of_all_observations 
             index_of_new_obs = partition_result$index_of_obs
             neighbor_names = partition_result$neighbor_names
             local_obs_list =  lapply(1:length(neighbor_names), function(i) voronoi_of_all_obs[neighbor_names[i]][[1]])
             if (index_of_new_obs %in% c(V_1_names[1], V_1_segments_names[1])){
               name_and_index = unlist(strsplit(index_of_new_obs, '_'))
               name_type = name_and_index[1]   ## = vertex
               ith = strtoi(name_and_index[2]) ## ith = 1
               initial_local_v = new_k_vertices[ith:(ith+1), ]
               update_v1_v2 = function(local_v, start_v = NULL, end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_v1_v2, initial_local_v, start_v = NULL, end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
             }else{
               ith = new_k_nbr
               initial_local_v = new_k_vertices[(ith-1):ith, ]
               update_2_vertices = function(local_v, start_v = NULL , end_v = NULL){
                 whole_vertices = data.matrix(rbind(start_v, local_v, end_v))
                 loss = 0
                 for (i in 1:length(neighbor_names)){
                   if (is.null(voronoi_of_all_obs[neighbor_names[i]][[1]])){
                     'pass'
                   }else{
                     
                     temp = name_index_seperation(neighbor_names[i])
                     if (temp[1] == 'segment'){
                       bg_index = strtoi(temp[2]) 
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2line(x, whole_vertices[bg_index,], whole_vertices[bg_index+1,])))
                       
                     }else{
                       bg_index = strtoi(temp[2])
                       local_loss = sum(apply(voronoi_of_all_obs[neighbor_names[i]][[1]], 1, function(x) squared_dist_point2point(x, whole_vertices[bg_index,])))
                       
                     }
                     loss = loss + local_loss
                   }
                 }
                 return(loss)
               }
               
               Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = NULL , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
               if (sum(is.na(Grd)) >0){
                 Grd = 0
               }
               
               new_local_v = initial_local_v - learning_rate * Grd  
               
               NN = 1
               while((max(new_local_v - initial_local_v) > epsilon)){
                 initial_local_v = new_local_v
                 Grd = matrix(grad(f = update_2_vertices, initial_local_v, start_v = NULL , end_v = NULL, "simple"), nrow = dim(initial_local_v)[1])
                 if (sum(is.na(Grd)) >0){
                   Grd = 0
                 }
                 
                 new_local_v = initial_local_v- learning_rate * Grd  
                 NN = NN+1
                 if(NN > 100){
                   break
                 }
                 
               }
               final_vertices = lgs(new_local_v, start_v = NULL, end_v = NULL, whole_lattice, n_nearest, inverse_eta_t)
               
             }
           },
           
           stop("Enter something that switches me!")
    )
    
    
    
    
    k_final_vertices[[as.character(k)]] = final_vertices  
    
  }
    
  vals = unlist(lapply(1:length(k_final_vertices), function(i) obt_function(k_final_vertices[[i]], rbind(past_obs, new_obs), inverse_eta_t, L, c1, c2, c3, p, R_local, delta)))
  
  
  return(k_final_vertices[[which(vals ==min(vals))[1]]])
  
  
}