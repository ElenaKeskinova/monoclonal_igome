rnd_graph = function(N,lett = 7,freqs=freq,s){
  
  a = sapply(1:N, function(x){
    paste(sample(sort(AA_STANDARD), lett, prob = freqs, replace = TRUE), collapse = "")
  })
  a = unique(a)
  g = newgraph(a,2)
  gl = graph_to_line(g)
  return(list(g,gl)) 
}