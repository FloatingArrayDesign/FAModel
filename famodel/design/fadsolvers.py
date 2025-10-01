# a file to hold the custom solvers used in FAD

import numpy as np
import matplotlib.pyplot as plt
import time
#from scipy.optimize import fsolve
#import scipy.optimize


# ================================ original above / modified below ===========================================


"""
def eval_func1(X, args):
    '''returns target outputs and also secondary outputs for constraint checks etc.'''

    # Step 1. break out design variables and arguments into nice names

    # Step 2. do the evaluation (this may change mutable things in args)

    # Step 3. group the outputs into objective function value and others

    return Y, oths



def step_func1(X, args, Y, oths, Ytarget, err, tol, iter, maxIter):
    '''General stepping functions, which can also contain special condition checks or other adjustments to the process

    '''

    # step 1. break out variables as needed

    # do stepping, as well as any conditional checks

    return dL   # returns dX (step to make)
"""



def dsolve1D(eval_func, step_func, X0, Ytarget, args, tol=0.0001, maxIter=20, Xmin=-np.inf, Xmax=np.inf):
    '''
        Assumes the function is positive sloped (so use -X if negative-sloped)

        tol        - relative convergence tolerance (relative to step size, dX)
        Xmin, Xmax - bounds. by default start bounds at infinity
    '''

    X = 1*X0         # start off design variable


    print(f"Starting dsolve1D iterations>>>   aiming for Y={Ytarget}")

    for iter in range(maxIter):


        # call evaluation function
        Y, oths = eval_func(X, args)

        # compute error
        err = Y - Ytarget

        print(f"  new iteration with X={X:6.2f} and Y={Y:6.2f}")

        # update/narrow the bounds (currently this part assumes that the function is positively sloped)  << any N-D equivalent?
        if   err > 0:# and L < LUpper:       #
            Xmax = 1.0*X
        elif err < 0:# and L > LLower:       #
            Xmin = 1.0*X

        if iter==maxIter-1:
            print("Failed to find solution after "+str(iter)+" iterations, with error of "+str(err))
            breakpoint()
            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>>
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?

        else:
            dX = step_func(X, args, Y, oths, Ytarget, err, tol, iter, maxIter)


        # check for convergence
        if np.abs(dX) < tol*(np.abs(X)+tol):
            print("Equilibrium solution completed after "+str(iter)+" iterations with error of "+str(err)+" and dX of "+str(dX))
            print("solution X is "+str(X))
            break


        # Make sure we're not diverging by keeping things within narrowing bounds that span the solution.
        #        I.e. detect potential for oscillation and avoid bouncing out and then back in to semi-taut config
        #        Use previous values to bound where the correct soln is, and if an iteration moves beyond that,
        #        stop it and put it between the last value and where the bound is (using golden ratio, why not).
        if   dX > 0 and X+dX >= Xmax:      # if moving up and about to go beyond previous too-high value
            X = X + 0.62*(Xmax-X)               # move to mid point between current value and previous too-high value, rather than overshooting
            print("<--|")
        elif dX < 0 and X+dX <= Xmin:      # if moving down and about to go beyond previous too-low value
            X = X + 0.62*(Xmin-X) #0.5*(L+LLower)               # move to mid point between current value and previous too-low value, rather than overshooting
            print("|-->")
        else:
            X = X+dX


    return X, Y, dict(iter=iter, err=err)



#    X, Y, info = dsolve1D(eval_func1, step_func1, X0, Ytarget, args, tol=tol, maxIter=maxIter)




# TODO: add default step_func (finite differencer), Ytarget, and args

def dsolve(eval_func, X0, Ytarget=[], step_func=None, args=[], tol=0.0001, maxIter=20,
           Xmin=[], Xmax=[], a_max=2.0, dX_last=[], display=0):
    '''
    PARAMETERS
    ----------
    eval_func : function
        function to solve (will be passed array X, and must return array Y of same size)
    X0 : array
        initial guess of X
    Ytarget : array (optional)
        target function results (Y), assumed zero if not provided
    stp_func : function (optional)
        function use for adjusting the variables (computing dX) each step.
        If not provided, Netwon's method with finite differencing is used.
    args : list
        A list of variables (e.g. the system object) to be passed to both the eval_func and step_func
    tol : float
        *relative* convergence tolerance (applied to step size components, dX)
    Xmin, Xmax
        Bounds. by default start bounds at infinity
    a_max
        maximum step size acceleration allowed
    dX_last
        Used if you want to dictate the initial step size/direction based on a previous attempt
    '''
    success = False

    # process inputs and format as arrays in case they aren't already

    X = np.array(X0, dtype=np.float_)         # start off design variable
    N = len(X)

    Xs = np.zeros([maxIter,N]) # make arrays to store X and error results of the solve
    Es = np.zeros([maxIter,N])
    dXlist = np.zeros([maxIter,N])
    dXlist2 = np.zeros([maxIter,N])


    # check the target Y value input
    if len(Ytarget)==N:
        Ytarget = np.array(Ytarget, dtype=np.float_)
    elif len(Ytarget)==0:
        Ytarget = np.zeros(N, dtype=np.float_)
    else:
        raise TypeError("Ytarget must be of same length as X0")


    # if a step function wasn't provided, provide a default one
    if step_func==None:
        if display>1:
            print("Using default finite difference step func")

        def step_func(X, args, Y, oths, Ytarget, err, tol, iter, maxIter):

            J = np.zeros([N,N])       # Initialize the Jacobian matrix that has to be a square matrix with nRows = len(X)

            for i in range(N):             # Newton's method: perturb each element of the X variable by a little, calculate the outputs from the
                X2 = np.array(X)                # minimizing function, find the difference and divide by the perturbation (finding dForce/d change in design variable)
                deltaX = tol*(np.abs(X[i])+tol)
                X2[i] += deltaX
                Y2, _, _ = eval_func(X2, args)    # here we use the provided eval_func

                J[:,i] = (Y2-Y)/deltaX             # and append that column to each respective column of the Jacobian matrix

            if N > 1:
                dX = -np.matmul(np.linalg.inv(J), Y-Ytarget)   # Take this nth output from the minimizing function and divide it by the jacobian (derivative)
            else:

                dX = np.array([-(Y[0]-Ytarget[0])/J[0,0]])

                if display > 1:
                    print(f" step_func iter {iter} X={X[0]:9.2e}, error={Y[0]-Ytarget[0]:9.2e}, slope={J[0,0]:9.2e}, dX={dX[0]:9.2e}")

            return dX                              # returns dX (step to make)



    # handle bounds
    if len(Xmin)==0:
        Xmin = np.zeros(N)-np.inf
    elif len(Xmin)==N:
        Xmin = np.array(Xmin, dtype=np.float_)
    else:
        raise TypeError("Xmin must be of same length as X0")

    if len(Xmax)==0:
        Xmax = np.zeros(N)+np.inf
    elif len(Xmax)==N:
        Xmax = np.array(Xmax, dtype=np.float_)
    else:
        raise TypeError("Xmax must be of same length as X0")



    if len(dX_last)==0:
        dX_last = np.zeros(N)
    else:
        dX_last = np.array(dX_last, dtype=np.float_)

    if display>1:
        print(f"Starting dsolve iterations>>>   aiming for Y={Ytarget}")


    for iter in range(maxIter):


        # call evaluation function
        Y, oths, stop = eval_func(X, args)

        # compute error
        err = Y - Ytarget

        if display>1:
            print(f"  new iteration #{iter} with X={X} and Y={Y}")

        Xs[iter,:] = X
        Es[iter,:] = err

        # stop if commanded by objective function
        if stop:
            break


        if iter==maxIter-1:
            if display>0:
                print("Failed to find solution after "+str(iter)+" iterations, with error of "+str(err))
                breakpoint()
            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>>
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?

        else:
            dX = step_func(X, args, Y, oths, Ytarget, err, tol, iter, maxIter)


        #if display>2:
        #    breakpoint()

        # Make sure we're not diverging by keeping things from reversing too much.
        # Track the previous step (dX_last) and if the current step reverses too much, stop it part way.
        # Stop it at a plane part way between the current X value and the previous X value (using golden ratio, why not).

        # get the point along the previous step vector where we'll draw the bounding hyperplane (could be a line, plane, or more in higher dimensions)
        Xlim = X - 0.62*dX_last

        # the equation for the plane we don't want to recross is then sum(X*dX_last) = sum(Xlim*dX_last)
        if np.sum((X+dX)*dX_last) < np.sum(Xlim*dX_last):         # if we cross are going to cross it

            alpha = np.sum((Xlim-X)*dX_last)/np.sum(dX*dX_last)    # this is how much we need to scale down dX to land on it rather than cross it

            if display > 2:
                print("  limiting oscillation with alpha="+str(alpha))
                print(f"   dX_last was {dX_last}, dX was going to be {dX}, now it'll be {alpha*dX}")
                print(f"   dX_last was {dX_last/1000}, dX was going to be {dX/1000}, now it'll be {alpha*dX/1000}")

            dX = alpha*dX  # scale down dX

        # also avoid extreme accelerations in the same direction
        if np.linalg.norm(dX_last) > tol:                           # only worry about accelerations if the last step was non-negligible
            for i in range(N):

                if abs(dX_last[i]) < tol:                           # set the maximum permissible dx in each direction based an an acceleration limit
                    dX_max = a_max*10*tol*np.sign(dX[i])
                else:
                    dX_max = a_max*dX_last[i]

                if dX_max == 0.0:                                   # avoid a divide-by-zero case (if dX[i] was zero to start with)
                    dX[i] = 0.0
                else:
                    a_i = dX[i]/dX_max                              # calculate ratio of desired dx to max dx

                    if a_i > 1.0:

                        if display > 2:
                            print(f"    limiting acceleration ({1.0/a_i:6.4f}) for axis {i}")
                            print(f"     dX_last was {dX_last}, dX was going to be {dX}")

                        #dX = dX*a_max/a_i  # scale it down to the maximum value
                        dX[i] = dX[i]/a_i  # scale it down to the maximum value (treat each DOF individually)

                        if display > 2:
                            print(f"     now dX will be {dX}")

        dXlist[iter,:] = dX
        if iter==196:
            breakpoint()
        # enforce bounds
        for i in range(N):

            if X[i] + dX[i] < Xmin[i]:
                dX[i] = Xmin[i] - X[i]

            elif X[i] + dX[i] > Xmax[i]:
                dX[i] = Xmax[i] - X[i]

        dXlist2[iter,:] = dX
        # check for convergence
        if all(np.abs(dX) < tol*(np.abs(X)+tol)):

            if display>0:
                print(f"dsolve converged. iter={iter}, X={X}, error={err} and dX={dX}")

                #if abs(err) > 10:
                #    breakpoint()

            if any(X == Xmin) or any(X == Xmax):
                success = False
                print("Warning: dsolve ended on a bound.")
            else:
                success = True

            break

        dX_last = 1.0*dX # remember this current value


        X = X + dX


    return X, Y, dict(iter=iter, err=err, dX=dX_last, oths=oths, Xs=Xs, Es=Es, success=success, dXlist=dXlist, dXlist2=dXlist2)


def dsolve2(eval_func, X0, Ytarget=[], step_func=None, args=[], tol=0.0001, maxIter=20,
           Xmin=[], Xmax=[], a_max=2.0, dX_last=[], stepfac=4, display=0):
    '''
    PARAMETERS
    ----------
    eval_func : function
        function to solve (will be passed array X, and must return array Y of same size)
    X0 : array
        initial guess of X
    Ytarget : array (optional)
        target function results (Y), assumed zero if not provided
    stp_func : function (optional)
        function use for adjusting the variables (computing dX) each step.
        If not provided, Netwon's method with finite differencing is used.
    args : list
        A list of variables (e.g. the system object) to be passed to both the eval_func and step_func
    tol : float or array
        If scalar, the*relative* convergence tolerance (applied to step size components, dX).
        If an array, must be same size as X, and specifies an absolute convergence threshold for each variable.
    Xmin, Xmax
        Bounds. by default start bounds at infinity
    a_max
        maximum step size acceleration allowed
    dX_last
        Used if you want to dictate the initial step size/direction based on a previous attempt
    '''
    success = False
    start_time = time.time()
    # process inputs and format as arrays in case they aren't already

    X = np.array(X0, dtype=np.float_)         # start off design variable
    N = len(X)

    Xs = np.zeros([maxIter,N]) # make arrays to store X and error results of the solve
    Es = np.zeros([maxIter,N])
    dXlist = np.zeros([maxIter,N])
    dXlist2 = np.zeros([maxIter,N])


    # check the target Y value input
    if len(Ytarget)==N:
        Ytarget = np.array(Ytarget, dtype=np.float_)
    elif len(Ytarget)==0:
        Ytarget = np.zeros(N, dtype=np.float_)
    else:
        raise TypeError("Ytarget must be of same length as X0")

    # ensure all tolerances are positive
    if np.isscalar(tol) and tol <= 0.0:
        raise ValueError('tol value passed to dsovle2 must be positive')
    elif not np.isscalar(tol) and any([toli <= 0 for toli in tol]):
        raise ValueError('every tol entry passed to dsovle2 must be positive')


    # handle bounds
    if len(Xmin)==0:
        Xmin = np.zeros(N)-np.inf
    elif len(Xmin)==N:
        Xmin = np.array(Xmin, dtype=np.float_)
    else:
        raise TypeError("Xmin must be of same length as X0")

    if len(Xmax)==0:
        Xmax = np.zeros(N)+np.inf
    elif len(Xmax)==N:
        Xmax = np.array(Xmax, dtype=np.float_)
    else:
        raise TypeError("Xmax must be of same length as X0")


    # if a step function wasn't provided, provide a default one
    if step_func==None:
        if display>1:
            print("Using default finite difference step func")

        def step_func(X, args, Y, oths, Ytarget, err, tols, iter, maxIter):
            ''' this now assumes tols passed in is a vector'''
            J = np.zeros([N,N])       # Initialize the Jacobian matrix that has to be a square matrix with nRows = len(X)

            for i in range(N):             # Newton's method: perturb each element of the X variable by a little, calculate the outputs from the
                X2 = np.array(X)                # minimizing function, find the difference and divide by the perturbation (finding dForce/d change in design variable)
                deltaX = stepfac*tols[i]                  # note: this function uses the tols variable that is computed in dsolve based on the tol input
                X2[i] += deltaX
                Y2, _, _ = eval_func(X2, args)    # here we use the provided eval_func

                J[:,i] = (Y2-Y)/deltaX             # and append that column to each respective column of the Jacobian matrix

            if N > 1:
                dX = -np.matmul(np.linalg.inv(J), Y-Ytarget)   # Take this nth output from the minimizing function and divide it by the jacobian (derivative)
            else:
                # if the result of the eval_func did not change, increase the stepfac parameter by a factor of 10 and calculate the Jacobian again
                if J[0,0] == 0.0:
                    
                    stepfacb = stepfac*10

                    J = np.zeros([N,N])       # Initialize the Jacobian matrix that has to be a square matrix with nRows = len(X)
                    for i in range(N):             # Newton's method: perturb each element of the X variable by a little, calculate the outputs from the
                        X2b = np.array(X)                # minimizing function, find the difference and divide by the perturbation (finding dForce/d change in design variable)
                        deltaXb = stepfacb*tols[i]                  # note: this function uses the tols variable that is computed in dsolve based on the tol input
                        X2b[i] += deltaXb
                        Y2b, _, _ = eval_func(X2b, args)    # here we use the provided eval_func
                        J[:,i] = (Y2b-Y)/deltaXb             # and append that column to each respective column of the Jacobian matrix
                    
                    if J[0,0] == 0.0:       # if the Jacobian is still 0, maybe increase the stepfac again, but there might be a separate issue
                        #breakpoint()
                        raise ValueError('dsolve2 found a zero gradient - maybe a larger stepfac is needed.')

                # if the Jacobian is all good, then calculate the dX
                dX = np.array([-(Y[0]-Ytarget[0])/J[0,0]])

                if display > 1:
                    print(f" step_func iter {iter} X={X[0]:9.2e}, error={Y[0]-Ytarget[0]:9.2e}, slope={J[0,0]:9.2e}, dX={dX[0]:9.2e}")

            return dX                              # returns dX (step to make)


    if len(dX_last)==0:
        dX_last = np.zeros(N)
    else:
        dX_last = np.array(dX_last, dtype=np.float_)

    if display>0:
        print(f"Starting dsolve iterations>>>   aiming for Y={Ytarget}")


    for iter in range(maxIter):


        # call evaluation function
        Y, oths, stop = eval_func(X, args)

        # compute error
        err = Y - Ytarget

        if display>2:
            print(f"  new iteration #{iter} with X={X} and Y={Y}")

        Xs[iter,:] = X
        Es[iter,:] = err

        # stop if commanded by objective function
        if stop:
            break

        # handle tolerances input
        if np.isscalar(tol):
            tols = tol*(np.abs(X)+tol)
        else:
            tols = np.array(tol)

        # check maximum iteration
        if iter==maxIter-1:
            if display>0:
                print("Failed to find solution after "+str(iter)+" iterations, with error of "+str(err))

            # looks like things didn't converge, so if N=1 do a linear fit on the last 30% of points to estimate the soln
            if N==1:

                m,b = np.polyfit(Es[int(0.7*iter):iter,0], Xs[int(0.7*iter):iter,0], 1)
                X = np.array([b])
                Y = np.array([0.0])
                if display>1:
                    print(f"Using linear fit to estimate solution at X={b}")

            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>>
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?

        else:
            dX = step_func(X, args, Y, oths, Ytarget, err, tols, iter, maxIter)


        #if display>2:
        #    breakpoint()

        # Make sure we're not diverging by keeping things from reversing too much.
        # Track the previous step (dX_last) and if the current step reverses too much, stop it part way.
        # Stop it at a plane part way between the current X value and the previous X value (using golden ratio, why not).

        # get the point along the previous step vector where we'll draw the bounding hyperplane (could be a line, plane, or more in higher dimensions)
        Xlim = X - 0.62*dX_last

        # the equation for the plane we don't want to recross is then sum(X*dX_last) = sum(Xlim*dX_last)
        if np.sum((X+dX)*dX_last) < np.sum(Xlim*dX_last):         # if we cross are going to cross it

            alpha = np.sum((Xlim-X)*dX_last)/np.sum(dX*dX_last)    # this is how much we need to scale down dX to land on it rather than cross it

            if display > 2:
                print("  limiting oscillation with alpha="+str(alpha))
                print(f"   dX_last was {dX_last}, dX was going to be {dX}, now it'll be {alpha*dX}")
                print(f"   dX_last was {dX_last/1000}, dX was going to be {dX/1000}, now it'll be {alpha*dX/1000}")

            dX = alpha*dX  # scale down dX

        # also avoid extreme accelerations in the same direction
        for i in range(N):

            if abs(dX_last[i]) > tols[i]:                           # only worry about accelerations if the last step was non-negligible

                dX_max = a_max*dX_last[i]                           # set the maximum permissible dx in each direction based an an acceleration limit

                if dX_max == 0.0:                                   # avoid a divide-by-zero case (if dX[i] was zero to start with)
                    breakpoint()
                    dX[i] = 0.0
                else:
                    a_i = dX[i]/dX_max                              # calculate ratio of desired dx to max dx

                    if a_i > 1.0:

                        if display > 2:
                            print(f"    limiting acceleration ({1.0/a_i:6.4f}) for axis {i}")
                            print(f"     dX_last was {dX_last}, dX was going to be {dX}")

                        #dX = dX*a_max/a_i  # scale it down to the maximum value
                        dX[i] = dX[i]/a_i  # scale it down to the maximum value (treat each DOF individually)

                        if display > 2:
                            print(f"     now dX will be {dX}")

        dXlist[iter,:] = dX
        #if iter==196:
            #breakpoint()

        # enforce bounds
        for i in range(N):

            if X[i] + dX[i] < Xmin[i]:
                dX[i] = Xmin[i] - X[i]

            elif X[i] + dX[i] > Xmax[i]:
                dX[i] = Xmax[i] - X[i]

        dXlist2[iter,:] = dX
        # check for convergence
        if all(np.abs(dX) < tols):

            if display>0:
                print("Iteration converged after "+str(iter)+" iterations with error of "+str(err)+" and dX of "+str(dX))
                print("Solution X is "+str(X))

                #if abs(err) > 10:
                #    breakpoint()

                if display > 0:
                    print("Total run time: {:8.2f} seconds = {:8.2f} minutes".format((time.time() - start_time),((time.time() - start_time)/60)))


            if any(X == Xmin) or any(X == Xmax):
                success = False
                print("Warning: dsolve ended on a bound.")
            else:
                success = True

            break

        dX_last = 1.0*dX # remember this current value


        X = X + dX


    return X, Y, dict(iter=iter, err=err, dX=dX_last, oths=oths, Xs=Xs, Es=Es, success=success, dXlist=dXlist, dXlist2=dXlist2)


def dsolvePlot(info):
    '''Plots dsolve or dsolve solution process based on based dict of dsolve output data'''

    n = info['Xs'].shape[1]  # number of variables

    if n < 8:
        fig, ax = plt.subplots(2*n, 1, sharex=True)
        for i in range(n):
            ax[  i].plot(info['Xs'][:info['iter']+1,i])
            ax[n+i].plot(info['Es'][:info['iter']+1,i])
        ax[-1].set_xlabel("iteration")
    else:
        fig, ax = plt.subplots(n, 2, sharex=True)
        for i in range(n):
            ax[i,0].plot(info['Xs'][:info['iter']+1,i])
            ax[i,1].plot(info['Es'][:info['iter']+1,i])
        ax[-1,0].set_xlabel("iteration, X")
        ax[-1,1].set_xlabel("iteration, Error")
    plt.show()


def dopt(eval_func, X0, tol=0.0001, maxIter=20, Xmin=[], Xmax=[], a_max=1.2, dX_last=[], display=0, stepfac=10):
    '''
        Multi-direction Newton's method solver.

        tol        - *relative* convergence tolerance (applied to step size components, dX)
        Xmin, Xmax - bounds. by default start bounds at infinity
        a_max      - maximum step size acceleration allowed
        stepfac    - factor to increase step size to relative to tol*X0
    '''
    start_time = time.time()

    success = False
    lastConverged = False # flag for whether the previous iteration satisfied the convergence criterion

    # process inputs and format as arrays in case they aren't already
    if len(X0) == 0:
        raise ValueError("X0 cannot be empty")

    X = np.array(X0, dtype=np.float_)         # start off design variable (optimized)

    # do a test call to see what size the results are
    f, g, Xextra, Yextra, oths, stop = eval_func(X) #, XtLast, Ytarget, args)

    N      = len(X)               # number of design variables
    Nextra = len(Xextra)          # additional relevant variables calculated internally and passed out, for tracking
    m      = len(g)               # number of constraints

    Xs = np.zeros([maxIter, N + Nextra]) # make arrays to store X and error results of the solve
    Fs = np.zeros([maxIter])     # make arrays to store objective function values
    Gs = np.zeros([maxIter, m]) # make arrays to store constraint function values



    if len(Xmin)==0:
        Xmin = np.zeros(N)-np.inf
    elif len(Xmin)==N:
        Xmin = np.array(Xmin, dtype=np.float_)
    else:
        raise TypeError("Xmin must be of same length as X0")

    if len(Xmax)==0:
        Xmax = np.zeros(N)+np.inf
    elif len(Xmax)==N:
        Xmax = np.array(Xmax, dtype=np.float_)
    else:
        raise TypeError("Xmax must be of same length as X0")


    if len(dX_last)==N:
        dX_last = np.array(dX_last, dtype=np.float_)
    elif len(dX_last)==0:
        dX_last = np.zeros(N)
    else:
        raise ValueError("dX_last input must be of same size as design vector, if provided")
    #XtLast = 1.0*Xt0

    # set finite difference step size
    #dX_fd = 4.0 #0.5# 1.0*dX[i]  # this is gradient finite difference step size, not opto step size
    dX_fd = stepfac*X*tol    # set dX_fd as function of tolerance and initial values



    if display > 0:
        print("Starting dopt iterations>>>")

    for iter in range(maxIter):
        iter_start_time = time.time()

        # call evaluation function (returns objective val, constrain vals, tuned variables, tuning results)
        f, g, Xextra, Yextra, oths, stop = eval_func(X) #, XtLast, Ytarget, args)

        if display > 1: print("")
        if display > 0:

            if isinstance(Xextra, list):
                XextraDisp = Xextra
            else:
                XextraDisp = Xextra.tolist()

            print((" >> Iteration {:3d}: f={:8.2e} X="+"".join([" {:9.2f}"]*len(X))+"  Xe="+"".join([" {:9.2f}"]*len(Xextra))).format(*(
                                 [ iter , f ]             + X.tolist()                  + XextraDisp)   ))


        if display > 1: print(f"\n    Constraint values: {g}")

        Xs[iter,:] = np.hstack([X, Xextra])
        Fs[iter]   = f
        Gs[iter,:] = g


        # stop if commanded by objective function
        if stop:
            message = 'Received stop command from objective function'
            break

        # temporarily display output
        #print(np.hstack([X,Y]))


        if iter==maxIter-1:

            print("Failed to converge after "+str(iter)+" iterations")

            if any(X == Xmin) or any(X == Xmax) or any(g < 0.0):
                for i in range(N):
                    if X[i] == Xmin[i] : print(f" Warning: Design variable {i} ended on minimum bound {Xmin[i]}.")
                    if X[i] == Xmax[i] : print(f" Warning: Design variable {i} ended on maximum bound {Xmax[i]}.")

                for j in range(m):                                 # go through each constraint
                    if g[j] < 0:                                     # if a constraint will be violated
                        print(f" Warning: Constraint {j} was violated by {-g[j]}.")
            else:
                print(" No constraint or bound issues.")

            success = False
            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>>
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?

        else:   # this is where we get derivatives and then take a step

            #dX = step_func(X, args, Y, oths, Ytarget, err, tol, iter, maxIter)
            # hard coding a generic approach for now

            dX = np.zeros(N)  # optimization step size to take

            X2 = np.array(X, dtype=np.float_)

            Jf = np.zeros([N])
            Jg = np.zeros([N,m])
            Hf = np.zeros([N])   # this is just the diagonal of the Hessians
            Hg = np.zeros([N,m])

            for i in range(N):                  # loop through each variable

                # could do repetition to hone in when second derivative is large, but not going to for now
                # or if first derivative is zero (in which case take a larger step size)

                X2[i] += dX_fd[i]                                 # perturb +
                fp, gp, Xtp, Yp, othsp, stopp = eval_func(X2)
                X2[i] -= 2.0*dX_fd[i]                             # perturb -
                fm, gm, Xtm, Ym, othsm, stopm = eval_func(X2)
                X2[i] += dX_fd[i]                                 # restore to original

                # for objective function and constraints (note that g may be multidimensional),
                # fill in diagonal values of Jacobian and Hession (not using off-diagonals for now)
                Jf[i]   = (fp-fm) /(2*dX_fd[i])
                Jg[i,:] = (gp-gm) /(2*dX_fd[i])
                Hf[i]   = (fm-2.0*f+fp) /dX_fd[i]**2
                Hg[i,:] = (gm-2.0*g+gp) /dX_fd[i]**2

            #breakpoint()

            # If we're currently violating a constraint, fix it rather than worrying about the objective function
            # This step is when new gradients need to be calculated at the violating point
            # e.g. in cases where the constraint functions are flat when not violated
            if any(g < 0.0):

                if display > 3:
                    print("         CONSTRAINT HANDLING SECTION")
                    for i in range(len(Jg)):
                        print(f"           Jg[{i}]     = {np.round(Jg[i],5)}")
                        #print(("               Jg[{:3d}] = "+"".join([" {:6.2f}"]*m).format(*([i]+Jg[i].tolist()))))

                g0    = []
                gradg = []
                #sqg   = []

                # first get the gradient of each active constraint
                stepdir = np.zeros(N)                               # this is the direction we will step in

                for j in range(m):                                  # go through each constraint
                    if g[j] < 0:                                    # if a constraint will be violated
                        if np.sum(np.abs(Jg[:,j])) == 0.0:
                            print(f"dopt error, zero Jacobian for constraint {j}. g(X) may be flat or dX_fd may be too small")
                            stop=True   # set flag to exit iteration
                            message = f"Error, zero Jacobian for constraint {j}. g(X) may be flat or dX_fd may be too small"
                            break

                        g0.append(   g[j])                    # constraint value at the current location
                        gradg.append(Jg[:,j])                 # gradient for each active constraint <<< doesn't work so well
                        #sqg.append( np.sum(Jg[:,j]*Jg[:,j]))  # gradient dotted with itself (i.e. sum of squares)


                        # OG output for comparison
                        stepdir_i = 1.0*Jg[:,j]   # default is to assume we're moving in the same direction as the gradient since that's most efficient
                        for i in range(N):
                            if (X[i]==Xmin[i] and Jg[i,j]<0) or (X[i]==Xmax[i] and Jg[i,j]>0):   # but if any dimension is on its bound, and the gradient is to move in that direction
                                stepdir_i[i] = 0.0                                                # set its component to zero instead (other dimensions will now have to move farther)
                        alph = (0.0-g[j])/np.sum(Jg[:,j]*stepdir_i)  # for our selected step direction, find how far to move to get to zero
                        if np.sum(Jg[:,j]*stepdir_i) == 0.0:
                            print('NaN isue')

                        dXcon = stepdir_i*alph *1.1                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint

                        if display > 3:
                            print(f'         - Looking at g[{j}]')
                            print("           stepdir_i = "+"".join([" {:.5f}"]*len(stepdir_i)).format(*(stepdir_i.tolist())))
                            print("           alph      = ",alph)
                            print("           g0        = ",g0)
                            print("           gradg     = ",gradg)

                        if display > 1:
                            print(("           Con {:3d} OG correction"+"".join([" {:9.2f}"]*N)).format(*(  [j]+ dXcon.tolist()) ))


                # now zero any dimensions that are about to cross a bound (if we're already at the bound)
                for i in range(N):

                    for j in range(len(g0)):   # look through each active constraint (but apply zeroing to all active constraints for now)
                        if (X[i]==Xmin[i] and gradg[j][i]<0) or (X[i]==Xmax[i] and gradg[j][i]>0):   # but if any dimension is on its bound, and the gradient is to move in that direction
                            for k in range(len(g0)):
                                gradg[k][i] = 0.0                                                # set its component to zero instead (other dimensions will now have to move farther)
                                if display > 3: print('gradg',gradg)
                if display > 3: print('         - No bounds issues')
                sqg = [ np.sum(jac*jac) for jac in gradg]  # update the gradient dotted with itself (i.e. sum of squares)


                if display > 3: print('         - Find stepdir')
                # now sort out a combined step direction depending on the active constraints
                if len(g0) == 2 and np.sum(gradg[0]*gradg[1]) < 0 and N>1:  # if two active constraints in opposing directions
                    c1 = g0[0]/sqg[0] * ( np.sum(gradg[0]*gradg[1]) * gradg[1]/sqg[1] - gradg[0] )

                    c2 = g0[1]/sqg[1] * ( np.sum(gradg[0]*gradg[1]) * gradg[0]/sqg[0] - gradg[1] )
                    stepdir = c1 + c2
                    if display > 3: print(f'           A: c1={c1}, c2={c2}')

                else:  # all other cases - assume we're moving in the same direction as the gradient since that's most efficient
                    #c2 = [(-g0[j])/np.sum(gradg[j]*gradg[j])*gradg[j] for j in range(len(g0))]    # compute step directions that will zero each constraint

                    c = np.zeros([len(g0), N])
                    for j in range(len(g0)):                          # compute step directions that will zero each constraint
                        if np.sum(gradg[j]*gradg[j]) > 0:                    # just leave it as zero if any direction has a zero derivative
                            c[j,:] = -g0[j] / np.sum(gradg[j]*gradg[j]) * gradg[j]
                            if display > 3: print(f'           B: c={c}')
                        else:
                            if display > 0:
                                print(f'  dopt warning: zero gradient squared for active constraint {j} at iter={iter} and X={X}')

                    #stepdir=sum(c2)
                    stepdir = np.sum(c, axis=0)                                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint
                if display > 3: print('           stepdir = ',stepdir)


                if np.linalg.norm(stepdir)==0:
                    stop = True
                    break


                if display > 3: print('         - Find alpha')
                # now find how large the step needs to be to satisfy each active constraint
                alpha = 0.0                                         # this is the scalar that determines how far we will step in the direction
                for j in range(m):                                  # go through each constraint
                    if g[j] < 0:                                    # if a constraint will be violated
                        alpha_i = (0.0-g[j])/np.sum(Jg[:,j]*stepdir)# for this constraint, find how far to move along the step direction to get to zero

                        alpha = np.max([alpha, alpha_i])
                        if display > 3: print('           alpha_i =',alpha_i)
                        # if an acceleration limit will be applied in some dimension, it'd be nice to revise the direction and recompute <<<

                        #dXcon = stepdir*alpha *1.1                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint

                        #if display > 1:
                            #print(f"    Constraint {j:3d} active).")
                            #print(("    Con {:3d} correction:  "+"".join([" {:9.2f}"]*N)).format(*(  [j]+ dXcon.tolist()) ))
                        #if display > 2:
                        #    print(("                            J = "+"".join([" {:9.2e}"]*m)).format(*Jg[:,j].tolist() ))
                        #    print(("                            H = "+"".join([" {:9.2e}"]*m)).format(*Hg[:,j].tolist() ))

                dX = stepdir*alpha *1.1       # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin)


                if display > 1:
                    print(("    Total constraint step (dX)         :"+"".join([" {:9.2f}"]*N)).format(*dX.tolist()) )

                #if iter==4 or iter==5:
                #    breakpoint()

                # if the above fails, we could try backtracking along dX_last until the constriant is no longer violated...

                # at the end of this, the step will be a summation of the steps estimated to resolve each constraint - good idea?

            # otherwise make an optimization step
            else:
                if display > 3:  print("         OPTIMIZATION STEP SECTION")
                
                # figure out step size in each dimension
                dxType = ['none']*N
                for i in range(N):
                    if Hf[i] <= 0.1*abs(Jf[i])/np.linalg.norm(dX_last):      # if the hessian is very small or negative, just move a fixed step size
                        #dX[i] = -Jf[i]/np.linalg.norm(Jf) * np.abs(dX_last[i]) * a_max*0.9
                        dX[i] = -Jf[i]/np.linalg.norm(Jf) * np.linalg.norm(dX_last) * a_max
                        if display > 3: print(Jf[i], np.linalg.norm(Jf), np.linalg.norm(dX_last), dX_fd[i])
                        # but make sure the step size is larger than the convergence tolerance
                        if abs(dX[i]) <= tol*(np.abs(X[i])+tol):
                           dX[i] = np.sign(dX[i])*tol*(np.abs(X[i])+tol)*1.1

                        dxType[i] = 'fixed'
                    else:
                        dX[i] = -Jf[i]/Hf[i]

                        dxType[i] = 'hessian'

                    #dX[i] = -Jf[i]/np.linalg.norm(Jf) * np.linalg.norm(dX_last) * a_max # << trying a fixed step size approach (no hessian)

                if display > 1:
                    print(("    Minimization step,     dX = "+"".join([" {:9.2f}"]*N)).format(*dX.tolist() ))
                if display > 2:
                    print(("                      step type "+"".join([" {:9}"]*N)).format(*dxType ))
                if display > 2:
                    print(("                            J = "+"".join([" {:9.2f}"]*N)).format(*Jf.tolist() ))
                    print(("                            H = "+"".join([" {:9.2f}"]*N)).format(*Hf.tolist() ))
                #breakpoint()

                if any(np.isnan(dX)):
                    breakpoint()

                dX_min0 = np.array(dX)

                # respect bounds (handle each dimension individually)
                for i in range(N):
                    if X[i] + dX[i] < Xmin[i]:
                        dX[i] = Xmin[i] - X[i]
                    elif X[i] + dX[i] > Xmax[i]:
                        dX[i] = Xmax[i] - X[i]

                dX_minB = np.array(dX)

                # deal with potential constraint violations in making the step (based on existing gradients)
                # respect constraints approximately (ignore cross-couplings...for now)
                X2 = X + dX                                         # save jump before constraint correction
                for j in range(m):                                  # go through each constraint
                    g2j = g[j] + np.sum(Jg[:,j]*dX)                 # estimate constraint value after planned step
                    if g2j < 0:                                     # if the constraint will be violated

                        # option 1: assume we complete the step, then follow the constraint gradient up to resolve the constraint violation
                        alpha = -g2j / np.sum(Jg[:,j]*Jg[:,j])      # assuming we follow the gradient, finding how far to move to get to zero

                        if display > 2 and alpha > 2000:
                            breakpoint()

                        dX = dX + alpha*Jg[:,j]*1.05                     # step size is gradient times alpha (adding a little extra for margin)


                        # option 2: just stop short of where the constraint would be violated along the original dX path (inferior because it gets bogged up when against a constraint)
                        #alpha = -g[j] / np.sum(Jg[:,j]*dX)
                        #dX = alpha * dX * 0.95

                        if display > 1:
                            print(("     trimin step: (j={:2d},{:6.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                                       *( [j,   alpha] +          dX.tolist())))


                # this is how to stop the dX vector at the approximate constraint boundary (not good for navigation)
                #for j in len(g):                           # go through each constraint
                #    if g[j] + np.sum(Jg[:,j]*dX) < 0:        # if a constraint will be violated
                #        alpha = -g[j]/np.sum(Jg[:,j]*dX)     # find where the constraint boundary is (linear approximation)
                #        dX = dX*alpha                        # shrink the step size accordingly (to stop at edge of constraint)

        if stop:
            break


        # Make sure we're not diverging by keeping things from reversing too much.
        # Track the previous step (dX_last) and if the current step reverses too much, stop it part way.

        '''
        # Original approach: Stop it at a plane part way between the current X value and the previous X value (using golden ratio, why not).
        # This means scaling down the full vector (while preserving its direction). The downside is this limits all dimensions.
        # get the point along the previous step vector where we could draw a bounding hyperplane (could be a line, plane, or more in higher dimensions)
        Xlim = X - 0.62*dX_last
        # the equation for the plane we don't want to recross is then sum(X*dX_last) = sum(Xlim*dX_last)
        if np.sum((X+dX)*dX_last) < np.sum(Xlim*dX_last):         # if we cross are going to cross it
            alpha = np.sum((Xlim-X)*dX_last)/np.sum(dX*dX_last)    # this is how much we need to scale down dX to land on it rather than cross it
            dX = alpha*dX  # scale down dX
            if display > 1:
                print(("                          (alpha={:9.2e})  to "+"".join([" {:8.2e}"]*N)).format(
        '''
        # Revised approach: only scale down the directions that have reversed sign from the last step.
        for i in range(N):
            if np.sign(dX[i])==-np.sign(dX_last[i]) and abs(dX_last[i]) > tol:   # if this dimension is reversing direction

                ratio = np.abs(dX[i]/dX_last[i])       # if it's reversing by more than 62% of the last step in this dimension, limit it
                if ratio > 0.62:
                    dX[i] = dX[i]*0.62/ratio           # scale down dX so that it is just 62% of the dX_last

                    if display > 1:
                        print(("     oscil limit: (i={:2d},{:6.3f},{:7.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                                     *( [i, 0.62/ratio, ratio]      +   dX.tolist() )))


        # also avoid extreme accelerations in the same direction
        if np.linalg.norm(dX_last) > tol:                           # only worry about accelerations if the last step was non-negligible
            for i in range(N):

                # set the maximum permissible dx in each direction based an an acceleration limit
                if abs(dX_last[i]) < tol:
                    dX_max = a_max*10*tol*np.sign(dX[i])
                #if abs(dX_last[i]) < tol*(np.abs(X[i])+tol):
                #    dX_max = a_max*tol*(np.abs(X[i])+tol)*np.sign(dX[i])
                else:
                    dX_max = a_max*dX_last[i]
                #print('dX_max',dX_max, dX_last, tol, dX)
                if dX_max == 0.0:                                   # avoid a divide-by-zero case (if dX[i] was zero to start with)
                    dX[i] = 0.0
                else:
                    a_i = dX[i]/dX_max                                  # calculate ratio of desired dx to max dx
                    #print('a_i',a_i, i, dX[i])
                    if a_i > 1.0:

                        #dX = dX/a_i  # scale it down to the maximum value   <<< this needs to be in the conditional  if X[i] > Xmin[i] and X[i] < Xmax[i]:  # limit this direction if it exceeds the limit and if it's not on a bound (otherwise things will get stuck)
                        # NOTE: if this has problems with making the dX too small, triggering convergence, try the individual approach below <<<
                        dX[i] = dX[i]/a_i  # scale it down to the maximum value (treat each DOF individually)
                        #print(dX[i])
                        if display > 1:
                            print(("     accel limit: (i={:2d},{:6.3f},{:7.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                        *( [i,   1.0/a_i, a_i ] +                   dX.tolist())))



        # enforce bounds
        for i in range(N):
            if X[i] + dX[i] < Xmin[i]:
                dX[i] = Xmin[i] - X[i]
                if display > 3: print(f"           Minimum bounds adjustment for dX[{i}]")
            elif X[i] + dX[i] > Xmax[i]:
                dX[i] = Xmax[i] - X[i]
                if display > 3: print(f"           Maximum bounds adjustment for dX[{i}]")


        # check for convergence
        if all(np.abs(dX) < tol*(np.abs(X)+tol)):

            if lastConverged:     # only consider things converged if the last iteration also satisfied the convergence criterion

                if display>0:
                    print(f"Optimization converged after {iter} iterations with dX of {dX}")
                    print(f"Solution X is "+str(X))
                    print(f"Constraints are "+str(g))

                if any(X == Xmin) or any(X == Xmax):
                    if display>0:
                        for i in range(N):
                            if X[i] == Xmin[i] : print(f" Warning: Design variable {i} ended on minimum bound {Xmin[i]}.")
                            if X[i] == Xmax[i] : print(f" Warning: Design variable {i} ended on maximum bound {Xmax[i]}.")

                    success = True
                    message = "converged on one or more bounds"

                elif any(X == Xmin) or any(X == Xmax) or any(g < 0.0):
                    if display>0:
                        for j in range(m):                                 # go through each constraint
                            if g[j] < 0:                                     # if a constraint will be violated
                                print(f" Warning: Constraint {j} was violated by {-g[j]}.")

                    success = False
                    message = f"converged with one or more constraints violated (by max {-min(g):7.1e})"

                else:
                    success = True
                    message = "converged with no constraint violations or active bounds"
                break

            else:
                lastConverged = True    #  if this is the first time the convergence criterion has been met, note it and keep going
                message = "convergence criteria only met once (need twice in a row)"
        else:
            lastConverged = False
            message = "not converged"

        if display > 2:
            print(f" Convergence message: {message}")

        dX_last = 1.0*dX # remember this current value
        #XtLast = 1.0*Xt
        #if iter==3:
            #breakpoint()
        X = X + dX

        if display > 1:

            print(("    dopt iteration finished. dX=        "+"".join([" {:9.2f}"]*N)).format(*(dX.tolist())))
            print("     iteration run time: {:9.2f} seconds".format(time.time() - iter_start_time))


    if display > 2:
        print(f" Convergence message: {message}")

    if display > 0:
        print("    total run time: {:8.2f} seconds = {:8.2f} minutes".format((time.time() - start_time),((time.time() - start_time)/60)))

    return X, f, dict(iter=iter, dX=dX_last, oths=oths, Xs=Xs, Fs=Fs, Gs=Gs, Xextra=Xextra, g=g, Yextra=Yextra,
                      success=success, message=message)



def dopt2(eval_func, X0, tol=0.0001, maxIter=20, Xmin=[], Xmax=[], a_max=1.2, dX_last=[], display=0, stepfac=10, args=[]):
    '''
        Gradient descent solver with some line search capability

        tol        - *relative* convergence tolerance (applied to step size components, dX)
        Xmin, Xmax - bounds. by default start bounds at infinity
        a_max      - maximum step size acceleration allowed
        stepfac    - factor to increase step size to relative to tol*X0
    '''
    start_time = time.time()

    success = False
    lastConverged = False # flag for whether the previous iteration satisfied the convergence criterion

    # process inputs and format as arrays in case they aren't already
    if len(X0) == 0:
        raise ValueError("X0 cannot be empty")

    X = np.array(X0, dtype=np.float_)         # start off design variable (optimized)

    # do a test call to see what size the results are
    f, g, Xextra, Yextra, oths, stop = eval_func(X, args) #, XtLast, Ytarget, args)

    N      = len(X)               # number of design variables
    Nextra = len(Xextra)          # additional relevant variables calculated internally and passed out, for tracking
    m      = len(g)               # number of constraints

    Xs = np.zeros([maxIter, N + Nextra]) # make arrays to store X and error results of the solve
    Fs = np.zeros([maxIter])     # make arrays to store objective function values
    Gs = np.zeros([maxIter, m]) # make arrays to store constraint function values



    if len(Xmin)==0:
        Xmin = np.zeros(N)-np.inf
    elif len(Xmin)==N:
        Xmin = np.array(Xmin, dtype=np.float_)
    else:
        raise TypeError("Xmin must be of same length as X0")

    if len(Xmax)==0:
        Xmax = np.zeros(N)+np.inf
    elif len(Xmax)==N:
        Xmax = np.array(Xmax, dtype=np.float_)
    else:
        raise TypeError("Xmax must be of same length as X0")


    if len(dX_last)==N:
        dX_last = np.array(dX_last, dtype=np.float_)
    elif len(dX_last)==0:
        dX_last = np.zeros(N)
    else:
        raise ValueError("dX_last input must be of same size as design vector, if provided")
    #XtLast = 1.0*Xt0

    # set finite difference step size
    #dX_fd = 4.0 #0.5# 1.0*dX[i]  # this is gradient finite difference step size, not opto step size
    dX_fd = stepfac*X*tol    # set dX_fd as function of tolerance and initial values
    dX_fd0 = np.array(dX_fd)

    if display > 3: print(f" dX_fd is {dX_fd}")

    if display > 0:
        print("Starting dopt iterations>>>")
    
    for iter in range(maxIter):
        iter_start_time = time.time()
        
        Xsave = np.array(X)

        if any(X == Xmin):
            Xbadj = np.array(X)
            for ixmin in np.where(X==Xmin)[0]:
                Xbadj[ixmin] = X[ixmin]*(1+tol)     # badj = bound adjustment
        elif any(X == Xmax):
            Xbadj = np.array(X)
            for ixmax in np.where(X==Xmax)[0]:
                Xbadj[ixmax] = X[ixmax]*(1-tol)
        else:
            Xbadj = np.array(X)
        
        X = np.array(Xbadj)
        
        # call evaluation function (returns objective val, constrain vals, tuned variables, tuning results)
        f, g, Xextra, Yextra, oths, stop = eval_func(X, args) #, XtLast, Ytarget, args)
        
        if display > 1: print("")
        if display > 0:

            if isinstance(Xextra, list):
                XextraDisp = Xextra
            else:
                XextraDisp = Xextra.tolist()

            print((" >> Iteration {:3d}: f={:8.2e} X="+"".join([" {:9.2f}"]*len(X))+"  Xe="+"".join([" {:9.2f}"]*len(Xextra))).format(*(
                                 [ iter , f ]             + X.tolist()                  + XextraDisp)   ))


        if display > 1: 
            print(f"\n    Constraint values: {g}")
        elif display > 0 and any(g < 0.0):
            print(f"    Constraint values: {g}")

        Xs[iter,:] = np.hstack([X, Xextra])
        Fs[iter]   = f
        Gs[iter,:] = g


        # stop if commanded by objective function
        if stop:
            message = 'Received stop command from objective function'
            break

        # temporarily display output
        #print(np.hstack([X,Y]))


        if iter==maxIter-1:

            print("Failed to converge after "+str(iter)+" iterations")

            if any(X == Xmin) or any(X == Xmax) or any(g < 0.0):
                for i in range(N):
                    if X[i] == Xmin[i] : print(f" Warning: Design variable {i} ended on minimum bound {Xmin[i]}.")
                    if X[i] == Xmax[i] : print(f" Warning: Design variable {i} ended on maximum bound {Xmax[i]}.")

                for j in range(m):                                 # go through each constraint
                    if g[j] < 0:                                     # if a constraint will be violated
                        print(f" Warning: Constraint {j} was violated by {-g[j]}.")
            else:
                print(" No constraint or bound issues.")

            success = False
            break

        #>>>> COULD ALSO HAVE AN ITERATION RESTART FUNCTION? >>>
        #  that returns a restart boolean, as well as what values to use to restart things if true. How?

        else:   # this is where we get derivatives and then take a step

            #dX = step_func(X, args, Y, oths, Ytarget, err, tol, iter, maxIter)
            # hard coding a generic approach for now

            dX = np.zeros(N)  # optimization step size to take

            X2 = np.array(X, dtype=np.float_)

            Jf = np.zeros([N])
            Jg = np.zeros([N,m])
            Hf = np.zeros([N])   # this is just the diagonal of the Hessians
            Hg = np.zeros([N,m])

            for i in range(N):                  # loop through each variable

                # could do repetition to hone in when second derivative is large, but not going to for now
                # or if first derivative is zero (in which case take a larger step size)

                dX_fd = np.array(dX_fd0)                        # make a copy of the original dX_fd to store temporary values
                
                X2[i] += dX_fd0[i]                              # perturb + by original dX_fd0
                if X2[i] > Xmax[i]:                             # if the perturbed+ X2 value goes above the bounds
                    X2[i] = Xmax[i]                             # set the perturbed+ X2 value to the max bound
                    dX_fd[i] = Xmax[i] - X[i]                     # and set the temp dX_fdi value to how much that new perturbation is
                
                fp, gp, Xtp, Yp, othsp, stopp = eval_func(X2, args)   # evaluate at the proper X2 position
                
                X2[i] -= 2.0*dX_fd[i]                           # perturb - by updated dX_fd
                if X2[i] < Xmin[i]:                             # if the perturbed- X2 value goes under the bounds
                    X2[i] = Xmin[i]                             # set the perturbed- X2 value to the min bound
                    dX_fd[i] = X[i] - Xmin[i]                     # and set the temp dX_fd value to how much that new perturbation is
                fm, gm, Xtm, Ym, othsm, stopm = eval_func(X2, args)   # evaluate at the proper X2 position
                
                X2[i] += dX_fd[i]                               # restore to original
                
                # for objective function and constraints (note that g may be multidimensional),
                # fill in diagonal values of Jacobian and Hessian (not using off-diagonals for now)
                Jf[i]   = (fp-fm) /(2*dX_fd[i])
                Jg[i,:] = (gp-gm) /(2*dX_fd[i])
                Hf[i]   = (fm-2.0*f+fp) /dX_fd[i]**2
                Hg[i,:] = (gm-2.0*g+gp) /dX_fd[i]**2
                #if i==0: print(fp, fm, dX_fd[i], Jf[i])
                
            #breakpoint()

            # If we're currently violating a constraint, fix it rather than worrying about the objective function
            # This step is when new gradients need to be calculated at the violating point
            # e.g. in cases where the constraint functions are flat when not violated
            if any(g < 0.0):

                if display > 3:
                    print("    Constraint step")
                    for i in range(len(Jg)):
                        print(f"           Jg[{i}]     = {np.round(Jg[i],5)}")
                        #print(("               Jg[{:3d}] = "+"".join([" {:6.2f}"]*m).format(*([i]+Jg[i].tolist()))))

                g0    = []
                gradg = []
                #sqg   = []

                # first get the gradient of each active constraint
                stepdir = np.zeros(N)                               # this is the direction we will step in

                for j in range(m):                                  # go through each constraint
                    if g[j] < 0:                                    # if a constraint will be violated
                        if np.sum(np.abs(Jg[:,j])) == 0.0:
                            print(f"dopt error, zero Jacobian for constraint {j}. g(X) may be flat or dX_fd may be too small")
                            stop=True   # set flag to exit iteration
                            message = f"Error, zero Jacobian for constraint {j}. g(X) may be flat or dX_fd may be too small"
                            break

                        g0.append(   g[j])                    # constraint value at the current location
                        gradg.append(Jg[:,j])                 # gradient for each active constraint <<< doesn't work so well
                        #sqg.append( np.sum(Jg[:,j]*Jg[:,j]))  # gradient dotted with itself (i.e. sum of squares)


                        # OG output for comparison
                        stepdir_i = 1.0*Jg[:,j]   # default is to assume we're moving in the same direction as the gradient since that's most efficient
                        for i in range(N):
                            if (X[i]==Xmin[i] and Jg[i,j]<0) or (X[i]==Xmax[i] and Jg[i,j]>0):   # but if any dimension is on its bound, and the gradient is to move in that direction
                                stepdir_i[i] = 0.0                                                # set its component to zero instead (other dimensions will now have to move farther)
                        alph = (0.0-g[j])/np.sum(Jg[:,j]*stepdir_i)  # for our selected step direction, find how far to move to get to zero
                        if np.sum(Jg[:,j]*stepdir_i) == 0.0:
                            print('NaN isue')

                        dXcon = stepdir_i*alph *1.1                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint

                        if display > 3:
                            print(f'         - Looking at g[{j}]')
                            print("           stepdir_i = "+"".join([" {:.5f}"]*len(stepdir_i)).format(*(stepdir_i.tolist())))
                            print("           alph      = ",alph)
                            print("           g0        = ",g0)
                            print("           gradg     = ",gradg)

                        if display > 1:
                            print(("           Con {:3d} OG correction"+"".join([" {:9.2f}"]*N)).format(*(  [j]+ dXcon.tolist()) ))


                # now zero any dimensions that are about to cross a bound (if we're already at the bound)
                for i in range(N):

                    for j in range(len(g0)):   # look through each active constraint (but apply zeroing to all active constraints for now)
                        if (X[i]==Xmin[i] and gradg[j][i]<0) or (X[i]==Xmax[i] and gradg[j][i]>0):   # but if any dimension is on its bound, and the gradient is to move in that direction
                            for k in range(len(g0)):
                                gradg[k][i] = 0.0                                                # set its component to zero instead (other dimensions will now have to move farther)
                                if display > 3: print('gradg',gradg)
                if display > 3: print('         - No bounds issues')
                sqg = [ np.sum(jac*jac) for jac in gradg]  # update the gradient dotted with itself (i.e. sum of squares)


                if display > 3: print('         - Find stepdir')
                # now sort out a combined step direction depending on the active constraints
                if len(g0) == 2 and np.sum(gradg[0]*gradg[1]) < 0 and N>1:  # if two active constraints in opposing directions
                    c1 = g0[0]/sqg[0] * ( np.sum(gradg[0]*gradg[1]) * gradg[1]/sqg[1] - gradg[0] )

                    c2 = g0[1]/sqg[1] * ( np.sum(gradg[0]*gradg[1]) * gradg[0]/sqg[0] - gradg[1] )
                    stepdir = c1 + c2
                    if display > 3: print(f'           A: c1={c1}, c2={c2}')

                else:  # all other cases - assume we're moving in the same direction as the gradient since that's most efficient
                    #c2 = [(-g0[j])/np.sum(gradg[j]*gradg[j])*gradg[j] for j in range(len(g0))]    # compute step directions that will zero each constraint

                    c = np.zeros([len(g0), N])
                    for j in range(len(g0)):                          # compute step directions that will zero each constraint
                        if np.sum(gradg[j]*gradg[j]) > 0:                    # just leave it as zero if any direction has a zero derivative
                            c[j,:] = -g0[j] / np.sum(gradg[j]*gradg[j]) * gradg[j]
                            if display > 3: print(f'           B: c={c}')
                        else:
                            if display > 0:
                                print(f'  dopt warning: zero gradient squared for active constraint {j} at iter={iter} and X={X}')

                    #stepdir=sum(c2)
                    stepdir = np.sum(c, axis=0)                                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint
                if display > 3: print('           stepdir = ',stepdir)


                if np.linalg.norm(stepdir)==0:
                    stop = True
                    break


                if display > 3: print('         - Find alpha')
                # now find how large the step needs to be to satisfy each active constraint
                alpha = 0.0                                         # this is the scalar that determines how far we will step in the direction
                for j in range(m):                                  # go through each constraint
                    if g[j] < 0:                                    # if a constraint will be violated
                        alpha_i = (0.0-g[j])/np.sum(Jg[:,j]*stepdir)# for this constraint, find how far to move along the step direction to get to zero

                        alpha = np.max([alpha, alpha_i])
                        if display > 3: print('           alpha_i =',alpha_i)
                        # if an acceleration limit will be applied in some dimension, it'd be nice to revise the direction and recompute <<<

                        #dXcon = stepdir*alpha *1.1                  # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin) - add the step command from each violated constraint

                        #if display > 1:
                            #print(f"    Constraint {j:3d} active).")
                            #print(("    Con {:3d} correction:  "+"".join([" {:9.2f}"]*N)).format(*(  [j]+ dXcon.tolist()) ))
                        #if display > 2:
                        #    print(("                            J = "+"".join([" {:9.2e}"]*m)).format(*Jg[:,j].tolist() ))
                        #    print(("                            H = "+"".join([" {:9.2e}"]*m)).format(*Hg[:,j].tolist() ))

                dX = stepdir*alpha *1.1       # step is step direction vector (possibly gradient) times alpha (plus a little extra for margin)


                if display > 1:
                    print(("    Total constraint step (dX)         :"+"".join([" {:9.2f}"]*N)).format(*dX.tolist()) )

                #if iter==4 or iter==5:
                #    breakpoint()

                # if the above fails, we could try backtracking along dX_last until the constriant is no longer violated...

                # at the end of this, the step will be a summation of the steps estimated to resolve each constraint - good idea?

            # otherwise (no constraints violated) make an optimization step
            else:
                
                # start by doing a line search down the slope
                
                dir = -Jf/np.linalg.norm(Jf)  # direction to move along                
                if display > 1: print(f"  beginning line search in steepest descent direction, u={dir}")
                step = dir * tol*np.linalg.norm(X) * 2
                #print('dir',dir,'step',step)
                j_active = -1                                       # index of constraint that is limiting the step
                                
                dX = np.zeros_like(X)
                X2 = X + dX 
                flast = 1.0*f
                glast = 1.0*g
                step2 = 1.0*step
                for k in range(100):  # now do a line search
                    
                    step2 = step*(2**k)   # double step size each time

                    dX = dX + step2   # <<< looks like I'm actually trippling the step - need to fix at some point <<<


                    # check for bound violation
                    if any(X2 + dX < Xmin) or any(X2 + dX > Xmax):  
                        dX = dX - step2                    
                        if display > 3:
                            print(f" ----- next step will cross bounds, so will use X={X+dX} and dX={dX}")
                        break

                    # evaluate the function
                    fl, gl, Xtl, Yl, othsl, stopl = eval_func(X + dX, args)   
                    if display > 3:
                        print(("      line searching: k={:2d} f={:6.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                               *( [k,   fl] +          (X+dX).tolist())))
                    
                    # check for increasing f
                    if fl > flast:
                        dX = dX - step2
                        if display > 3:
                            print(f" ----- f increasing ----- so will use X={X+dX} and dX={dX}")
                        break
                    
                    # check for constraint violation
                    if any(gl < 0):
                        
                        frac = -glast/(gl-glast)           # how much of the step (fraction) could be taken until each constraint is violated
                        backfrac = (1.0-frac)*(frac > 0)   # how much of the step to backtrace (with filtering to exclude negatives from constraints that are decreasing)               
                        j_active = np.argmax(backfrac)
                        
                        dXog = 1.0*dX
                        #breakpoint()
                        
                        dX = dX - backfrac[j_active]*step2 - 2*dir*tol*np.linalg.norm(X)  # back up, and also keep a 2*tol margin from the boundary
                        
                        # normal case: # find -Jf component tangent with constraint surface to move along
                        tandir = -Jf + np.sum(Jg[:,j_active]*Jf)*Jg[:,j_active]/np.sum(Jg[:,j_active]*Jg[:,j_active]) 
                    
                        '''   ...not sure this next part really works/helps at all...
                        # >>>>>>> make it so that if we've backed up further from the constraint than the last X, 
                        #then move along direction of previous dX! (to avoid potential sawtooth bouncing along constraint boundary)
                        # IF the constraint direction points more toward the constraint than the previous dX direction.
                        if iter > 0 and np.sum(dX*step) < 0:  # line search has us moving further away from constraint boundary
                            
                            tandir = tandir/np.linalg.norm(tandir)
                            lastdir = dX_last/np.linalg.norm(dX_last)
                            
                            if np.sum(Jg[:,j_active]*lastdir) > np.sum(Jg[:,j_active]*tandir):
                                print(f"special case. Using {lastdir} rather than {tandir}")
                                tandir = lastdir
                        '''    
                        
                        
                        if display > 3:
                            print(f" ----- constraint violated ----- {gl} ")
                            print(f"    will back up to X={X+dX} and do line search along constraint")
                                 
                            print(dXog)
                            print(dX)
                            print(gl)
                        fl2, gl2, Xtl2, Yl2, othsl2, stopl2 = eval_func(X + dX, args)
                        if display > 3: print(gl2)
                        
                        break
                        
                    
                    flast = 1.0*fl
                    glast = 1.0*gl
                
                #for i in range(N):
                    #dX[i] = -Jf[i]/np.linalg.norm(Jf) * np.linalg.norm(dX_last) * a_max # << trying a fixed step size approach (no hessian)

                if display > 1:
                    print(("    Minimization step,     dX = "+"".join([" {:9.2f}"]*N)).format(*dX.tolist() ))
                if display > 2:
                    print(("                            J = "+"".join([" {:9.2f}"]*N)).format(*Jf.tolist() ))

                if any(np.isnan(dX)):
                    breakpoint()

                dX_min0 = np.array(dX)

                # respect bounds (handle each dimension individually) <<<
                for i in range(N):
                    if X[i] + dX[i] < Xmin[i]:
                        dX[i] = Xmin[i] - X[i]
                    elif X[i] + dX[i] > Xmax[i]:
                        dX[i] = Xmax[i] - X[i]

                dX_minB = np.array(dX)

                
                # but do a line search tangent to whatever constraint boundary is limiting if applicable (and if tangent is clear)
                if j_active >= 0 and N > 1 and not any(np.isnan(tandir)):
                    #tandir = -Jf + np.sum(Jg[:,j_active]*Jf)*Jg[:,j_active]/np.sum(Jg[:,j_active]*Jg[:,j_active]) # find -Jf component tangent with constraint surface to move along
                    step = tandir/np.linalg.norm(tandir) * tol*np.linalg.norm(X) 
                    if display > 3:
                        print(f"Constraint normal vector is {Jg[:,j_active]/np.linalg.norm(Jg[:,j_active])}")
                        print(f"  beginning line search along constraint {j_active} boundary, u={tandir/np.linalg.norm(tandir)}")
                    
                    X3 = X + dX 
                    step3=0
                    for k in range(100):  # now do a line search
                        
                        # evaluate the function
                        fl, gl, Xtl, Yl, othsl, stopl = eval_func(X3, args)   
                        if display > 3:
                            print(("      line searching: k={:2d} f={:6.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                                   *( [k,   fl] +          X3.tolist())))
                        
                        # check for increasing f
                        if k>0:
                            if fl > flast:
                                X3 = X3 - step3
                                if display > 3:
                                    print(f" ----- f increasing ----- so will use previous X={X3} and dX={X3 - X}")
                                break
                                
                        # check for constraint violation
                        if any(gl < 0):
                            X3 = X3 - step3
                            # could instead back up to an intermediate point, and offset by the 2*tol margin too
                            if display > 3:
                                print(f" ----- constraint violated ----- {gl} --- so will use previous")
                            break
                            
                        flast = 1.0*fl
                        step3 = step*(1.6**k)  # increase step size each time
                        
                        # check for bound violation
                        if any(X3 + step3 < Xmin) or any(X3 + step3 > Xmax):                        
                            if display > 3:
                                print(f" ----- next step will cross bounds, so stopping here")
                            break
                        
                        X3 = X3 + step3   
                    
                    dX = X3 - X   # undo the last step (which was bad) and calculated overall effective dX                        


                # this is how to stop the dX vector at the approximate constraint boundary (not good for navigation)
                #for j in len(g):                           # go through each constraint
                #    if g[j] + np.sum(Jg[:,j]*dX) < 0:        # if a constraint will be violated
                #        alpha = -g[j]/np.sum(Jg[:,j]*dX)     # find where the constraint boundary is (linear approximation)
                #        dX = dX*alpha                        # shrink the step size accordingly (to stop at edge of constraint)

        if stop:
            break


        # Make sure we're not diverging by keeping things from reversing too much.
        # Track the previous step (dX_last) and if the current step reverses too much, stop it part way.

        
        # Original approach: Stop it at a plane part way between the current X value and the previous X value (using golden ratio, why not).
        # This means scaling down the full vector (while preserving its direction). The downside is this limits all dimensions.
        # get the point along the previous step vector where we could draw a bounding hyperplane (could be a line, plane, or more in higher dimensions)
        Xlim = X - 0.62*dX_last
        # the equation for the plane we don't want to recross is then sum(X*dX_last) = sum(Xlim*dX_last)
        if np.sum((X+dX)*dX_last) < np.sum(Xlim*dX_last):         # if we cross are going to cross it
            ratio = np.sum((Xlim-X)*dX_last)/np.sum(dX*dX_last)    # this is how much we need to scale down dX to land on it rather than cross it
            dX = ratio*dX  # scale down dX
            if display > 1:
                print(("     oscil limit: ( reducing by factor {:6.3f} :"+"".join([" {:9.2f}"]*N)).format(
                                                     *( [ratio]      +   dX.tolist() )))

        # Revised approach: only scale down the directions that have reversed sign from the last step.
        '''
        for i in range(N):
            if np.sign(dX[i])==-np.sign(dX_last[i]) and abs(dX_last[i]) > tol:   # if this dimension is reversing direction

                ratio = np.abs(dX[i]/dX_last[i])       # if it's reversing by more than 62% of the last step in this dimension, limit it
                if ratio > 0.62:
                    dX[i] = dX[i]*0.62/ratio           # scale down dX so that it is just 62% of the dX_last

                    if display > 1:
                        print(("     oscil limit: (i={:2d},{:6.3f},{:7.3f}):"+"".join([" {:9.2f}"]*N)).format(
                                                     *( [i, 0.62/ratio, ratio]      +   dX.tolist() )))
        '''

        # also avoid extreme accelerations in the same direction
        if np.linalg.norm(dX_last) > tol:                           # only worry about accelerations if the last step was non-negligible
            for i in range(N):

                # set the maximum permissible dx in each direction based an an acceleration limit
                #if abs(dX_last[i]) < tol:
                #    dX_max = a_max*10*tol*np.sign(dX[i])
                if abs(dX_last[i]) < tol*(np.abs(X[i])+tol):
                    dX_max = a_max*tol*(np.abs(X[i])+tol)*np.sign(dX[i])
                else:
                    dX_max = a_max*dX_last[i]
                #print('dX_max',dX_max, dX_last, tol, dX)
                if dX_max == 0.0:                                   # avoid a divide-by-zero case (if dX[i] was zero to start with)
                    dX[i] = 0.0
                else:
                    a_i = dX[i]/dX_max                                  # calculate ratio of desired dx to max dx
                    #print('a_i',a_i, i, dX[i])
                    if a_i > 1.0:

                        # Option 1. the directoin-preserving approach: (could have problems with making the dX too small, triggering convergence)
                        dX = dX/a_i  # scale it down to the maximum value   <<< this needs to be in the conditional  if X[i] > Xmin[i] and X[i] < Xmax[i]:  # limit this direction if it exceeds the limit and if it's not on a bound (otherwise things will get stuck)
                        # Option 2. the individual approach below <<<
                        #dX[i] = dX[i]/a_i  # scale it down to the maximum value (treat each DOF individually)
                        #print(dX[i])
                        if display > 1:
                            print(("     accel limit: (i={:2d},by {:6.3f}   :"+"".join([" {:9.2f}"]*N)).format(
                                        *( [i,   1.0/a_i] +                   dX.tolist())))
                        

        # enforce bounds
        for i in range(N):
            if X[i] + dX[i] < Xmin[i]:
                dX[i] = Xmin[i] - X[i]
                #dX[i] = Xmin[i]*(1+tol) - X[i]
                if display > 2: print(f"           Minimum bounds adjustment for dX[{i}]")
            elif X[i] + dX[i] > Xmax[i]:
                dX[i] = Xmax[i] - X[i]
                #dX[i] = Xmax[i]*(1-tol) - X[i]
                if display > 2: print(f"           Maximum bounds adjustment for dX[{i}]")


        # check for convergence
        if all(np.abs(dX) < tol*(np.abs(X)+tol)):

            if lastConverged:     # only consider things converged if the last iteration also satisfied the convergence criterion

                if display>0:
                    print(f"Optimization converged after {iter} iterations with dX of {dX}")
                    print(f"Solution X is "+str(X))
                    print(f"Constraints are "+str(g))

                if any(X == Xmin) or any(X == Xmax):
                    if display>0:
                        for i in range(N):
                            if X[i] == Xmin[i] : print(f" Warning: Design variable {i} ended on minimum bound {Xmin[i]}.")
                            if X[i] == Xmax[i] : print(f" Warning: Design variable {i} ended on maximum bound {Xmax[i]}.")

                    success = True
                    message = "converged on one or more bounds"

                elif any(X == Xmin) or any(X == Xmax) or any(g < 0.0):
                    if display>0:
                        for j in range(m):                                 # go through each constraint
                            if g[j] < 0:                                     # if a constraint will be violated
                                print(f" Warning: Constraint {j} was violated by {-g[j]}.")

                    success = False
                    message = f"converged with one or more constraints violated (by max {-min(g):7.1e})"

                else:
                    success = True
                    message = "converged with no constraint violations or active bounds"
                break

            else:
                lastConverged = True    #  if this is the first time the convergence criterion has been met, note it and keep going
                message = "convergence criteria only met once (need twice in a row)"
        else:
            lastConverged = False
            message = "not converged"

        if display > 2:
            print(f" Convergence message: {message}")

        dX_last = 1.0*dX # remember this current value
        #XtLast = 1.0*Xt
        #if iter==3:
            #breakpoint()
        X = X + dX

        if display > 0:
            print(("    dopt iteration finished. dX= "+"".join([" {:9.2f}"]*N)).format(*(dX.tolist())))
        if display > 2:
            print("     iteration run time: {:9.2f} seconds".format(time.time() - iter_start_time))


    if display > 2:
        print(f" Convergence message: {message}")

    if display > 0:
        print("    total run time: {:8.2f} seconds = {:8.2f} minutes".format((time.time() - start_time),((time.time() - start_time)/60)))
    
    runtime = time.time() - start_time  #seconds
    
    return X, f, dict(iter=iter, dX=dX_last, oths=oths, Xs=Xs, Fs=Fs, Gs=Gs, Xextra=Xextra, g=g, Yextra=Yextra,
                      success=success, message=message, time=runtime)


def doptPlot(info):

    n = info['Xs'].shape[1] # number of DVs
    m = info['Gs'].shape[1] # number of constraints

    fig, ax = plt.subplots(n+1+m,1, sharex=True)
    Xs = np.array(info["Xs"])
    Fs = np.array(info["Fs"])
    Gs = np.array(info["Gs"])
    iter = info["iter"]

    for i in range(n):
        ax[i].plot(Xs[:iter+1,i])

    ax[n].plot(Fs[:iter+1])
    ax[n].set_ylabel("objective")

    for i in range(Gs.shape[1]):
        j = i+1+n
        ax[j].axhline(0, color=[0.5,0.5,0.5])
        ax[j].plot(Gs[:iter+1,i])
        ax[j].set_ylabel(f'con {i}')

    ax[j].set_xlabel("iteration")

    plt.show()



# ------------------------------ sample functions ----------------------------


def eval_func1(X, args):
    '''returns target outputs and also secondary outputs for constraint checks etc.'''

    # Step 1. break out design variables and arguments into nice names

    # Step 2. do the evaluation (this may change mutable things in args)
    y1 = (X[0]-2)**2 + X[1]
    y2 = X[0] + X[1]

    # Step 3. group the outputs into objective function value and others
    Y = np.array([y1, y2])               # objective function
    oths = dict(status=1)                # other outputs - returned as dict for easy use

    return Y, oths, False



def step_func1(X, args, Y, oths, Ytarget, err, tol, iter, maxIter):
    '''General stepping functions, which can also contain special condition checks or other adjustments to the process

    '''

    # get numerical derivative
    J = np.zeros([len(X),len(X)])       # Initialize the Jacobian matrix that has to be a square matrix with nRows = len(X)

    for i in range(len(X)):             # Newton's method: perturb each element of the X variable by a little, calculate the outputs from the
        X2 = np.array(X)                # minimizing function, find the difference and divide by the perturbation (finding dForce/d change in design variable)
        deltaX = tol*(np.abs(X[i])+tol)
        X2[i] += deltaX
        Y2, extra = eval_func1(X2, args)

        J[:,i] = (Y2-Y)/deltaX             # and append that column to each respective column of the Jacobian matrix

    dX = -np.matmul(np.linalg.inv(J), Y)   # Take this nth output from the minimizing function and divide it by the jacobian (derivative)

    return dX                              # returns dX (step to make)





## ============================== below is a new attempt at the catenary solve ======================================
# <<< moved to Catenary.py >>>












'''

# test run


#Catenary2(100, 50, 130, 1e8, 100, plots=1)

print("\nTEST 1")
catenary(576.2346666666667, 514.6666666666666, 800, 4809884.623076923, -2.6132152062554828, CB=-64.33333333333337, HF0=0, VF0=0, Tol=1e-05, MaxIter=50, plots=2)
print("\nTEST 2")
catenary(88.91360441490338, 44.99537159734132, 100.0, 854000000.0000001, 1707.0544275185273, CB=0.0, HF0=912082.6820817506, VF0=603513.100376363, Tol=1e-06, MaxIter=50, plots=1)
print("\nTEST 3")
catenary(99.81149090002897, 0.8459770263789324, 100.0, 854000000.0000001, 1707.0544275185273, CB=0.0, HF0=323638.97834178555, VF0=30602.023233123222, Tol=1e-06, MaxIter=50, plots=1)
print("\nTEST 4")
catenary(99.81520776134033, 0.872357398602503, 100.0, 854000000.0000001, 1707.0544275185273, CB=0.0, HF0=355255.0943810993, VF0=32555.18285808794, Tol=1e-06, MaxIter=50, plots=1)
print("\nTEST 5")
catenary(99.81149195956499, 0.8459747131565791, 100.0, 854000000.0000001, 1707.0544275185273, CB=0.0, HF0=323645.55876751675, VF0=30602.27072107738, Tol=1e-06, MaxIter=50, plots=1)
print("\nTEST 6")
catenary(88.91360650151807, 44.99537139684605, 100.0, 854000000.0000001, 1707.0544275185273, CB=0.0, HF0=912082.6820817146, VF0=603513.100376342, Tol=1e-06, MaxIter=50, plots=1)
'''
'''
maxIter = 10
# call the master solver function
X0      = [2,2]
Ytarget = [0,0]
args    = []
X, Y, info = dsolve(eval_func1, step_func1, X0, Ytarget, args, maxIter=maxIter)
'''
