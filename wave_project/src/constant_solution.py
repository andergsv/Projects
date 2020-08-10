from solvers import *

def constant_solution(c):

    f = lambda x,y,t: 0
    q = lambda x,y: 2.
    I = lambda x,y: c
    V = lambda x,y: 0
    ue = lambda x, y, t: c
    
    diff =scalar(I, V, f, q, 2, 10, 50, 10, 50, 0.01, 1, 
                            ue, user_action='diff2', plot=False) 
    return diff


def test_constant():
    c = 0.8
    #ue = lambda x, y, t: c
    diff = constant_solution(c)
    tol = 1e-12
    assert diff < tol
    print 'Success'
   
test_constant()
