rnd_graph = function(N,lett = 7,freqs=freq,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  a = unique(a)
  g = newgraph(a,2)
  gl = graph_to_line(g)
  return(list(g,gl)) 
}



rnd_graph_reps = function(N,lett = 7,freqs=freq,repfreqs,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  a = unique(a)
  reps = sample(as.numeric(names(repfreqs)),length(a),prob = repfreqs,replace = T)
  names(reps) = a
  g = newgraph(a,2)
  g = set_vertex_attr(g,name = "Freq",value = reps)
  gl = graph_to_line(g)
  w = weights(gl,gl$edg_lcs,reps)
  gl = set_edge_attr(gl,name = "weight",value = w)
  
  return(list(g,gl)) 
}