"""A few simple math functions"""
# todo : numpy might be faster

def median(X):
    """Calculate median value given list of values."""
    X.sort()
    n = len(X)
    mid = (n - 1) / 2.0
    if n == 0:
        return None
    elif (n - 1) % 2 != 0:
        try:
            m = (X[int(mid)] + X[int(mid + 1)]) / 2.0
        except IndexError:
            print n
            m = None
        except TypeError:
            print mid, X[int(mid)], X[int(mid + 1)], len(X)
    else:
        m = float(X[int(mid)])
    return m
    
def mean(X):
    """Calculate mean of list X."""
    try:
        x = float(sum(X)) / float(len(X))
    except ZeroDivisionError:
        x = None
    return x

def variance(X, sample=False):
    """Calculate variance of list X."""
    try:
        m = float(sum(X)) / float(len(X))
    except ZeroDivisionError:
        return 0
    v_sum = 0.0
    for x in X:
        v_sum += (x - m)**2
    if not sample:
        v = v_sum / len(X)
    else:
        v = v_sum / (len(X) - 1)
    
    return v

def stdev(X, sample=False):
    """Calculate standard deviation of vector X."""
    s = variance(X, sample=sample) ** 0.5
    return s
    
def percentile(X, p):
    """Calculate percentile value given list of values."""
    X.sort()
    N = len(X)
    if N == 0:
        return None
    elif p < 1.0:
        n = p * (N - 0.5)
        try:
            m = float(X[int(round(n))])
        except TypeError:
            print 'Warning: None value returned for percentile'
            m = None
    else:
        m = float(X[-1])
    return m
    
def correlation(X, Y, sample=False):
    """Pearson correlation coefficient between vectors X & Y."""
    if len(X) != len(Y) or len(X) == 0:
        return None
    else:
        n = len(X)
        mx = mean(X)
        my = mean(Y)
        sx = stdev(X, sample=sample)
        sy = stdev(Y, sample=sample)
        sum = 0.0
        if sx * sy != 0:
            for i in xrange(0, n):
                sum += ((X[i] - mx) / sx) * ((Y[i] - my) / sy)
            if not sample:
                r = (1.0 / n) * sum
            else:
                r = (1.0 / (n - 1.0)) * sum
        else:
            return None
    return r