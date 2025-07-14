rnd_graph = function(N,lett = 7,freqs=freq,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  a = unique(a)
  g = newgraph(a,2)
  gl = graph_to_line(g)
  return(list(g,gl)) 
}


rnd_graph_w = function(N,lett = 7,freqs=freq,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  fr = table(a)
  f = fr[which((fr)>1)]
  aa = names(fr[which((fr)>1)])
  g = newgraph(aa,2)
  g = set_vertex_attr(g,name="Freq",value = f)
  gl = graph_to_line(g)
  
  weights = weights(gl,gl$edg_lcs,f)
  gl = set_edge_attr(gl,name = "weight",value = weights)
  
  return(list(g,gl)) 
}

rnd_graph_wcopy = function(N,lett = 7,freqs=freq,pepfr,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  
  g = newgraph(a,2)
  g = set_vertex_attr(g,name="Freq",value =pepfr)
  gl = graph_to_line(g)
  
  names(pepfr) = V(g)$name
  pepfr = log2(pepfr)
  weights = weights(gl,gl$edg_lcs,pepfr)
  gl = set_edge_attr(gl,name = "weight",value = weights)
  
  return(list(g,gl)) 
}
