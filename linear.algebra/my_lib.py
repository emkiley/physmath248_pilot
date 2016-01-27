def dot_prod(u,v):
    '''
    vector product
    
    input: two vectors of equal lenght
    output: scalar 
        
    '''
    if len(u) is not len(v): 
        print "error"
    uv_prod = 0
    for i in range(len(u)):
        uv_prod += u[i]*v[i]
    return uv_prod

