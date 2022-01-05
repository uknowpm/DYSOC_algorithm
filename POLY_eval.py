import numpy as np

def POLY_eval(X, b, flag):
    '''uses polynomial regression model to predict objective function values
    --------------------------------------------------------------------------

    Input 
    S - matrix containing points where prediction wanted
    b - parameter vector computed with POLY.m
    flag - order of polynomial. must be same as flag  used for computing vector b

    Output
    Yest - predicted objective function values corresponding to the points in
    S
    --------------------------------------------------------------------------
    '''
    if len(X.shape)==1:
        m=1
        n=X.shape[0]
        ones = np.ones(m)
        axis = 0
    else:
        m, n = X.shape
        ones = np.ones((m,1))
        axis = 1
    
    # linear
    if flag == 'lin':
        Xs = np.concatenate((ones, X),axis)
        Yest = np.dot(Xs, b)
    elif flag == 'quad':
        Xs = np.concatenate((ones, np.concatenate((X, np.multiply(X,X)),axis)),axis)
        #Xs = np.concatenate((np.ones((m, 1)), np.concatenate((X, X**2), axis = 1)), axis = 1)
        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                if m!=1:
                    x = np.multiply(X[:, ii].ravel() , X[:, j].ravel())#S[:, ii] * S[:, j]
                else:
                    if len(X.shape) == 1:
                        x = np.multiply(X[ii], X[j])
                    else:
                        x = np.multiply(X[0][ii] , X[0][j])
                #x = X[:, ii] * X[:, j]
                j = j + 1
                if len(X.shape) != 1:
                    Xs = np.concatenate((Xs, x.T),axis)
                else:
                    Xs = np.concatenate((Xs,np.array([x])),axis)
            ii = ii + 1
        Yest = np.dot(Xs, b)
    elif flag == 'cub':
        X_temp = np.concatenate((ones, np.concatenate((X, np.multiply(X,X)),axis)),axis)
        #X_temp = np.concatenate((np.ones((m, 1)), np.concatenate((X, X**2), axis = 1)), axis = 1)
        Xs = np.concatenate((X_temp, np.multiply(X,np.multiply(X,X))),axis)
        #Xs = np.concatenate((X_temp, X**3), axis = 1)
        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                if m!=1:
                    x = np.multiply(X[:, ii].ravel() , X[:, j].ravel())#S[:, ii] * S[:, j]
                else:
                    if len(X.shape) == 1:
                        x = np.multiply(X[ii] , X[j])
                    else:
                        x = np.multiply(X[0][ii], X[0][j])
                #x = np.array([np.multiply(X[:, ii] , X[:, j])])
                #x = X[:, ii] * X[:, j]
                j = j + 1
               
                if len(X.shape) != 1:
                    Xs = np.concatenate((Xs, np.matrix(x).T),axis)
                else:
                    Xs = np.concatenate((Xs,np.array([x])),axis)
                #Xs = np.concatenate((Xs, x), axis = 0)
            ii = ii + 1

        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                kk = j + 1
                while kk <= n - 1:
                    if m!=1:
                        x = np.multiply(X[:, ii].ravel() , X[:, j].ravel())#S[:, ii] * S[:, j]
                        #x = np.multiply(X[:, ii].ravel() , np.multiply(X[:, j].ravel(),X[:,kk].ravel()))
                    else:
                        if len(X.shape) == 1:
                            x = np.multiply(X[ii],X[j])
                        else:
                            x = np.multiply(X[0][ii] , X[0][j])
                        #x = np.multiply(X[0][ii] , np.multiply(X[0][j],X[0][kk]))
                        
                    #x = X[:, ii] * X[:, j] * X[:, kk]
                    kk = kk + 1
                    if len(X.shape) !=1:
                        Xs = np.concatenate((Xs, np.matrix(x).T),axis)
                    else:
                        Xs = np.concatenate((Xs,np.array([x])),axis)
                    #Xs = np.concatenate((Xs, x), axis = 1)
                j = j + 1
            ii = ii + 1
        Yest = np.dot(Xs, b)
    elif flag == 'quadr':
        if X.shape[0]==X.size:
            ones=ones[0]
        Xs = np.concatenate((ones, np.concatenate((X, np.multiply(X,X)),axis=1)),axis)
        #X = np.concatenate((np.ones((m, 1)), np.concatenate((X, X**2), axis = 1)), axis = 1)
        Yest = np.dot(Xs, b)
    elif flag == 'cubr':
        if X.shape[0]==X.size:
            ones=ones[0]
        X_temp = np.concatenate((ones, np.concatenate((X, np.multiply(X,X)),axis=1)),axis)
        Xs = np.concatenate((X_temp, np.multiply(X,np.multiply(X,X))),axis)
        #X_temp = np.concatenate((np.ones((m, 1)), np.concatenate((X, X**2), axis = 1)), axis = 1)
        #Xs = np.concatenate((X_temp, X**3), axis = 1)
        Yest = np.dot(Xs, b)

    Yest_ret=Yest.T
    Yest_ret=np.asarray(Yest_ret)
       
    if len(Yest_ret)==1 or Yest_ret.shape[0]==1:
         Yest_ret=np.array(np.matrix(Yest_ret).T)

    return Yest_ret