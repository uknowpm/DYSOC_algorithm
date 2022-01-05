import numpy as np

def POLY(S, Y, flag):
    '''computes the coefficients of a regression polynomial of the specified type
    --------------------------------------------------------------------------
    Input:
    S - sample site matrix
    Y - objective function values corresponding to the points in S
    flag - determines the order of polynomial model

    Output:
    b - vector with parameters
    --------------------------------------------------------------------------
    '''
    m, n = S.shape
    
    # linear
    if flag == 'lin':
        X = np.concatenate((np.ones((m, 1)), S), axis = 1)
        temp = X.conj().transpose()
        b = np.dot(np.linalg.solve(np.dot(temp, X), X.conj().transpose()), Y)
    elif flag == 'quad':
        X = np.concatenate((np.ones((m, 1)), np.concatenate((S, np.multiply(S,S)), axis = 1)), axis = 1)
        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                x = np.array([np.multiply(S[:, ii] , S[:, j])])#S[:, ii] * S[:, j]
                j = j + 1
                #X = np.concatenate((X, np.asmatrix(x).T), axis = 1)
                X = np.concatenate((X, x.T), axis = 1)
            ii = ii + 1
        temp = X.conj().transpose()
        b = np.dot(np.linalg.solve(np.dot(temp, X), X.conj().transpose()), Y)
    elif flag == 'cub':
       
        X_temp = np.concatenate((np.ones((m, 1)), np.concatenate((S, np.multiply(S,S)), axis = 1)), axis = 1)
        X = np.concatenate((X_temp, np.multiply(S,np.multiply(S,S))), axis = 1)
        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                x = np.array([np.multiply(S[:, ii] , S[:, j])])
                j = j + 1
                X = np.concatenate((X, x.T), axis = 1)
            ii = ii + 1

        ii = 0
        while ii < n - 1:
            j = ii + 1
            while j <= n - 1:
                kk = j + 1
                while kk <= n - 1:
                    
                    x = np.array([np.multiply(S[:, ii] , np.multiply(S[:, j],S[:,kk]))])
                    #x = S[:, ii] * S[:, j] * S[:, kk]
                    kk = kk + 1
                    X = np.concatenate((X, x.T), axis = 1)
                j = j + 1
            ii = ii + 1
        temp = X.conj().transpose()
     
        b = np.dot(np.linalg.solve(np.dot(temp, X), temp), Y)
     
        
        #print 'linear system result',np.dot( np.dot(X.T, X),np.linalg.solve(np.dot(temp, X), temp))
        print 'b',b.shape
    
        
    elif flag == 'quadr':
        X = np.concatenate((np.ones((m, 1)), np.concatenate((S, np.multiply(S,S)), axis = 1)), axis = 1)
        temp = X.conj().transpose()
        b = np.dot(np.linalg.solve(np.dot(temp, X), X.conj().transpose()), Y)
    elif flag == 'cubr':
        X_temp = np.concatenate((np.ones((m, 1)), np.concatenate((S, np.multiply(S,S)), axis = 1)), axis = 1)
        X = np.concatenate((X_temp, np.multiply(S,np.multiply(S,S))), axis = 1)
        temp = X.conj().transpose()
        b = np.dot(np.linalg.solve(np.dot(temp, X), X.conj().transpose()), Y)

    return b