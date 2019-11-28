### basic principal curve functions

## calculate squared distance between a point and line segement formed by v1 and v2
squared_dist_point2line = function(obs, v1, v2, angle_ind = FALSE){
  temp = sum((obs-v1)*(v2-v1))
  if (angle_ind){
    if (temp > 0 & sum((obs-v2)*(v1-v2))>0){ ## case when the projection of obs to line lies on the line
      distance = sum((obs-v1)^2) - (temp/sqrt(sum((v2-v1)^2)))^2  
    }else{distance = 10^6} #set this value big enough since we care about the minimum distance of an obs to Voronoi-like partition in function get_closest_Voronoi_index_and_distance
  }else{
  distance = sum((obs-v1)^2) - (temp/sqrt(sum((v2-v1)^2)))^2    
  }
  return (distance)
}




## calculate squared distance between obs point and vertice point
squared_dist_point2point = function(obs, v){
  return (sum((obs-v)^2)) 
}




## calculate Voronoi-like partition given a polygonal line and an obs, return partition index and distance of this obs to polygonal line
get_closest_Voronoi_index_and_distance = function(vertices, obs){
  nbr_vertices = length(vertices[,1]) ## number of vertices of polygonal lines. 
  dist2vertices = apply(vertices, 1, function(x) squared_dist_point2point(obs, x))
  dist2segments = apply(matrix(1:(nbr_vertices-1), ncol=1), 1, function(x) squared_dist_point2line(obs, vertices[x,], vertices[x+1,], TRUE))
  closest_vertices_index = which(dist2vertices == min(dist2vertices))[1]
  if(nbr_vertices == 1){
    print("The number of vertices should be strictly larger than 1 !")
  }
  if(nbr_vertices >2){
    if (closest_vertices_index == 1 & sum((obs-vertices[1,])*(vertices[2,] - vertices[1,])) <= 0){
      index = 'vertex_1'
      distance_of_obs = dist2vertices[1]
    }else if((closest_vertices_index == nbr_vertices) & (sum((obs-vertices[nbr_vertices,])*(vertices[(nbr_vertices-1),] - vertices[nbr_vertices,])) <= 0)){
      index = paste('vertex_', nbr_vertices, sep='')
      distance_of_obs = dist2vertices[nbr_vertices]
    }else if(closest_vertices_index %in% 2:(nbr_vertices-1)){  
      if ((sum((obs-vertices[closest_vertices_index,])*(vertices[(closest_vertices_index-1),] - vertices[closest_vertices_index,])) <= 0) & (sum((obs-vertices[closest_vertices_index,])*(vertices[(closest_vertices_index +1),] - vertices[closest_vertices_index,])) <= 0)){
          index = paste('vertex_', closest_vertices_index, sep='')
          distance_of_obs = dist2vertices[closest_vertices_index]
          }else{
            closest_segments_index = which(dist2segments == min(dist2segments))[1]
            if (min(dist2vertices) > min(dist2segments)){
              index = paste('segment_', closest_segments_index, sep='')
              distance_of_obs = dist2segments[closest_segments_index]
            }else{
              index = paste('vertex_', closest_vertices_index, sep='') 
              distance_of_obs = dist2vertices[closest_vertices_index]
            }
          }
    }else{
      closest_segments_index = which(dist2segments == min(dist2segments))[1]
      if (min(dist2vertices) > min(dist2segments)){
        index = paste('segment_', closest_segments_index, sep='')
        distance_of_obs = dist2segments[closest_segments_index]
      }else{
        index = paste('vertex_', closest_vertices_index, sep='')
        distance_of_obs = dist2vertices[closest_vertices_index]
      }
    }
  }else{
    if (closest_vertices_index == 1 & sum((obs-vertices[1,])*(vertices[2,] - vertices[1,])) <= 0){
      index = 'vertex_1'
      distance_of_obs = dist2vertices[1]
    }else if((closest_vertices_index == nbr_vertices) & (sum((obs-vertices[nbr_vertices,])*(vertices[(nbr_vertices-1),] - vertices[nbr_vertices,])) <= 0)){
      index = paste('vertex_', nbr_vertices, sep='')
      distance_of_obs = dist2vertices[nbr_vertices]
    }else{
      closest_segments_index = which(dist2segments == min(dist2segments))[1]
      if (min(dist2vertices) > min(dist2segments)){
        index = paste('segment_', closest_segments_index, sep='')
        distance_of_obs = dist2segments[closest_segments_index]
      }else{
        index = paste('vertex_', closest_vertices_index, sep='') 
        distance_of_obs = dist2vertices[closest_vertices_index]
      }
    } 
  }
  return (list('index'= index, 'distance_of_obs' = distance_of_obs))  
}




## given the vertices, find closest partition index for each obs in data set past_obs, and the neighbors of closest partition index of new_obs
voronoi_partitions = function(vertices, past_obs, new_obs){
  d = length(new_obs)
  vertices = matrix(vertices, ncol = d)
  nbr_vertices = length(vertices[,1])
  vertices_names = sapply(1:nbr_vertices, function(x) paste('vertex_', x, sep=''))
  segments_names = sapply(1:(nbr_vertices-1), function(x) paste('segment_', x, sep=''))
  #vertices_names = apply(matrix(1:nbr_vertices, ncol=1), 1, function(x) paste('vertex_', x, sep=''))
  #segments_names = apply(matrix(1:(nbr_vertices-1), ncol=1), 1, function(x) paste('segment_', x, sep=''))
  dictionary = vector('list', 2*nbr_vertices-1)  ## creat a list with vertices and segments name and obs                                                      index attached to them.
  names(dictionary) = c(vertices_names, segments_names)
  
  for (i in 1:length(past_obs[,1])){
    obs = past_obs[i,]
    id = get_closest_Voronoi_index_and_distance(vertices, obs)[[1]] ##index id = vertice_i or segment_i
    dictionary[id][[1]] = rbind(dictionary[id][[1]], obs)
  }
  
  index_of_new_obs = get_closest_Voronoi_index_and_distance(vertices, new_obs)[[1]]
  dictionary[index_of_new_obs][[1]] = rbind(dictionary[index_of_new_obs][[1]], new_obs)
  
  if (index_of_new_obs == vertices_names[1]){
    neighbor_names = c(vertices_names[1:2], segments_names[1:2])
  }else if (index_of_new_obs == vertices_names[nbr_vertices]){
    neighbor_names = c(vertices_names[(nbr_vertices-1):nbr_vertices], segments_names[(nbr_vertices-2):(nbr_vertices-1)])
  }else if (index_of_new_obs %in% segments_names){
    i_th = which(segments_names == index_of_new_obs)
    if (i_th == length(segments_names)){
      neighbor_names = c(vertices_names[i_th:(i_th+1)], segments_names[(i_th-1):i_th])
    }else{neighbor_names = c(vertices_names[i_th:(i_th+1)], segments_names[(i_th-1):(i_th+1)])}
    
  }else{
    i_th = which(vertices_names == index_of_new_obs)
    neighbor_names = c(vertices_names[(i_th-1):(i_th+1)], segments_names[(i_th-2):(i_th+1)])
  }

  return (list('partition_of_all_observations' = dictionary, 'index_of_obs' = index_of_new_obs, 'neighbor_names' = sort(neighbor_names)))
}




## add, keep or delete a vertex when one has a new observation, and return an updated vertices
add_keep_delete_vertex = function(vertices, new_k, new_obs){
  #vertices = matrix(vertices, ncol = length(new_obs))
  nbr_vertices = length(vertices[,1])
  vertices_names = sapply(1:nbr_vertices, function(x) paste('vertex_', x, sep=''))
  segments_names = sapply(1:(nbr_vertices-1), function(x) paste('segment_', x, sep=''))
  #vertices_names = unlist(lapply(1:nbr_vertices, function(x) paste('vertex_', x, sep='')))
  #segments_names = unlist(lapply(1:(nbr_vertices-1), function(x) paste('segment_', x, sep='')))
  
  if (new_k == (nbr_vertices+1)){
    decision = 'add'
  }else if (new_k == nbr_vertices){
    decision = 'keep'
  }else if(new_k == (nbr_vertices -1)){
    decision = 'delete'
  }else{
    decision = 'error'
  }
  index_of_new_obs = get_closest_Voronoi_index_and_distance(vertices, new_obs)[[1]]
  switch(decision,
    add = {
      if (index_of_new_obs == vertices_names[1]){
        new_vertex = new_obs
        new_vertices = rbind(new_obs, vertices)
        #neighbors_names = c(vertices_names[1], segments_names[1]) 
      }else if (index_of_new_obs == vertices_names[nbr_vertices]){
        new_vertex = new_obs
        new_vertices = rbind(vertices, new_obs) 
        #neighbors_names = c(vertices_names[nbr_vertices], segments_names[length(segments_names)])
      }else if (index_of_new_obs %in% segments_names){
        i_th = which(segments_names == index_of_new_obs)
        new_vertex = (vertices[i_th,] + vertices[(i_th+1),])/2
        new_vertices = rbind(vertices[1:i_th, ], new_vertex, vertices[(i_th+1):nbr_vertices,])
        #neighbors_names = c(vertices_names[i_th:(i_th+1)], segments_names[(i_th-1):(i_th+1)])
      }else{
        i_th = which(vertices_names == index_of_new_obs)
        new_vertex = (vertices[i_th,] + vertices[(i_th+1),])/2
        new_vertices = rbind(vertices[1:i_th, ], new_vertex, vertices[(i_th+1):nbr_vertices,])
        #neighbors_names = c(vertices_names[(i_th-1):(i_th+1)], segments_names[(i_th-2):(i_th+1)])
      }
    },
    
    keep = {
      new_vertices = vertices
    },
    
    delete = {
      if (index_of_new_obs == segments_names[1]){
        new_vertices = vertices[-2,]
      }else if(index_of_new_obs == vertices_names[1]){
        new_vertices = rbind(new_obs, vertices[-c(1,2), ])
      }else if (index_of_new_obs == vertices_names[nbr_vertices]){
        new_vertices = rbind(vertices[1:(nbr_vertices-1),], new_obs)   
      }else if(index_of_new_obs == segments_names[nbr_vertices-1]){
        new_vertices = vertices[-(nbr_vertices-1), ]
      }else if (index_of_new_obs %in% segments_names){
        i_th = which(segments_names == index_of_new_obs)
        new_vertices = vertices[-(i_th+1), ]
      }else{
        i_th = which(vertices_names == index_of_new_obs)
        new_vertices = vertices[-i_th, ]
      }
    },
    stop("Enter something that switches me!")
    )
  return (data.matrix(new_vertices)) 
}


## name index seperation function, will be used in OLPC_gamma.R
name_index_seperation = function(x){
  temp = unlist(strsplit(x, '_'))
  return (temp)
}




## cumulative loss function, the third parameter lambda is useless
cumul_loss = function(all_obs, current_vertices, lambda=0){
  current_squared_distances = (apply(all_obs, 1, function(x) get_closest_Voronoi_index_and_distance(current_vertices, x)[[2]]))
  #past_squared_distances = unlist(lapply(1:length(all_obs[,1]), function(x) get_closest_Voronoi_index_and_distance(past_vertices[[x]], all_obs[x])[[2]]))
  biais = sum(current_squared_distances)
  #variance = lambda/2* sum((current_squared_distances - past_squared_distances)^2)
  return (biais)
}




## prior over different k
log_prior_k = function(k, R, L, dlt, d){
  value = (k+1)*d/2*log(pi) - (k+1)*log(gamma(d/2+1)) + d*log(2*sqrt(d)*R) - d*log(dlt) + k*d*(log(L)-log(k*dlt)) 
  return (10*value) ##????? multiplicator 10
}


## length of polygonal line
length_of_polygonal_line = function(vertices){
  if (dim(vertices)[1] == 2){
    return (sqrt(sum((vertices[2:(dim(vertices)[1]) ,] - vertices[1:(dim(vertices)[1]-1) ,])^2)))
  }else{
    return (sum(sqrt(rowSums((vertices[2:(dim(vertices)[1]) ,] - vertices[1:(dim(vertices)[1]-1) ,])^2))))
}
}


## generate lattice within a certain area. The side length of lattice is delta.
lattice_generator = function(local_data, delta){
  local_center = colMeans(local_data)
  semi_diameter  = sqrt(max(rowSums((local_data - local_center)^2)))
  d = dim(local_data)[2]
  lll = lapply(1:d, function(x) seq(min(local_data[, x]), max(local_data[, x]), by = min(delta, (max(local_data[, x])- min(local_data[, x]))/5))) #make sure the lattice could be dense enough
  lttc = expand.grid(lll)
  # seems range and 'dlt' is useless in the future ?????
  range = unlist(lapply(1:d, function(x) (max(local_data[, x])- min(local_data[, x]))/5))
  return (list("lattice" = lttc, 'diameter' = 2*semi_diameter/sqrt(d), 'dlt' = min(delta, min(range))))
}





## penalty function; indicator shows whether to use P, the maximum number of vertices allowed, in penalty function. 
h = function(vertices, L, c1, c2, c3, p, R, delta, indicator = FALSE){
  if (indicator){
    k = p
  }else{
  k = length(vertices[,1]) #number of vertices
  }
  d = length(vertices[1,])  #dimension of vertices
  l = length_of_polygonal_line(vertices) #length of polygonal line
  V_d = pi^(d/2)/gamma(d/2+1) 
  term_1 = c1*(log(8*p*exp(1)*V_d)+3*d^(1.5)-d)*k
  term_2 = c2*(log(2)/sqrt(d)/delta + d/delta) * l
  term_3 = c3*d*(log(sqrt(d)*(2*R+delta))-log(delta))
  return (term_1+term_2+term_3)
}



## whole lattice with range R
lattice_generator_R = function(R, delta, d){
  lll = lapply(1:d, function(x) seq(-R*sqrt(d), R*sqrt(d), by = delta)) 
  lttc = expand.grid(lll)
  return (data.matrix(lttc))
}



## generate lattice along each dimension
lattice_generator_dimensions = function(mydata, delta, d ){
  lll = lapply(1:d, function(x) seq(min(mydata[, x]), max(mydata[, x]), by = delta))
  lttc = expand.grid(lll)
  return(data.matrix(lttc))
  
}


## given a data point, find_k_nearest_vertices
find_k_nearest_vertices = function(lttc, query, nn = 3){
  idex = nn2(lttc, query, nn)$nn.idx
  idex = unique(as.vector(idex))
  return (lttc[idex, ])
}

## minimising function
#object_function = function(whole_vertices_temp, obs_available, inverse_eta_t, L, c1,c2,c3,p,R_local,delta){
#   whole_vertices_temp = na.omit(whole_vertices_temp)
  
#   if (dim(whole_vertices_temp)[1] > 1){
#  loss_temp = cumul_loss(obs_available, whole_vertices_temp) + inverse_eta_t*(h(whole_vertices_temp, L, c1, c2, c3, p, R_local, delta)+sample(c(-1,1), 1, prob = c(0.5,0.5))*rexp(1,1))}else{loss_temp = Inf}
#return (loss_temp)
#  }

### 

obt_function = function(whole_vertices_temp, obs_available, inverse_eta_t, L, c1, c2, c3, p, R_local, delta){
  whole_vertices_temp = data.matrix(na.omit(whole_vertices_temp))
  loss_temp_1 = cumul_loss(obs_available, whole_vertices_temp, whole_vertices_temp) 
  loss_temp_2 = inverse_eta_t*(h(whole_vertices_temp, L, c1, c2, c3, p, R_local, delta) + sample(c(-1,1), 1, prob = c(0.5,0.5))*rexp(1,1)) #h is the penalty function
  return (loss_temp_1+loss_temp_2)
}





### local greedy search

lgs = function(local_v, start_v = NULL, end_v = NULL, whole_lattice, n_nearest, inverse_eta_t){
  nbr = dim(local_v)[1]
  d = dim(local_v)[2]
  # given a set of local vertices(local_v), for each of them, find the most n_nearest points from a lattice, and compose them as a list.
  # i.e., locally update vertices.
  local_lattice_list = lapply(1:nbr, function(x) find_k_nearest_vertices(whole_lattice, matrix(local_v[x,], nrow=1), nn =n_nearest))
  
  # find all possible combinations of locally updated vertices(at most 3) and combine each of them with the start and the end of polygonal line vertice
  # to form a new polygonal line
  local_vertices_list = list()
  if (nbr == 3){ # when 3 vertices need to be updated
    count = 1
     for (x in 1:n_nearest){
       v1 = local_lattice_list[[1]][x,]
       for (y in 1:n_nearest){
        v2 = local_lattice_list[[2]][y,]
         for (z in 1:n_nearest){
           v3 = local_lattice_list[[3]][z,]
          local_vertices_list[[count]] = data.matrix(rbind(start_v, v1, v2, v3, end_v))
          count = count + 1
         }
       }
     }
    # compute the loss for each updated polygonal line and return updated vertices of polygonal lines that attains the minimum
    rrr = unlist(lapply(1:n_nearest^nbr, function(x) obt_function(local_vertices_list[[x]], rbind(past_obs, new_obs), inverse_eta_t, L, c1, c2, c3, p, R_local, delta)))
    return (local_vertices_list[[which(rrr==min(rrr))[1]]])
    
    }
  
  if (nbr == 2){ # when 2 vertices need to be updated
    count = 1
    for (x in 1:n_nearest){
      v1 = local_lattice_list[[1]][x,]
      for (y in 1:n_nearest){
        v2 = local_lattice_list[[2]][y,]
        local_vertices_list[[count]] = data.matrix(rbind(start_v, v1, v2, end_v))
        count = count + 1
        
      }
    }
    rrr = unlist(lapply(1:n_nearest^nbr, function(x) obt_function(local_vertices_list[[x]], rbind(past_obs, new_obs), inverse_eta_t, L, c1, c2, c3, p, R_local, delta)))
    return (local_vertices_list[[which(rrr==min(rrr))[1]]])
  }
}
  
  
  
  
  
  
  
  
 







