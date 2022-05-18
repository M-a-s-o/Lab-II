import sys
import numpy as np

sys.path.append('../Scripts/utils/')
import propagate as pr

def main():
    print("Calibrazione")
    x=pr.GenerateXs(2)
    y=pr.GenerateXs(1)
    rappr = x[0]*632.8e-9/(2*x[1]*1e-6)
    rapprmean = [28, 10]
    rapprcov = np.diag([(1/np.sqrt(3))**2]*2)
    rapprres = pr.Propagate(rapprmean, rapprcov, 2, rappr, x)
    print(f"Rapporto: {rapprres}")
    prod = x[0]*x[1]
    prodmean = [rapprres[0], 10]
    prodcov = np.diag([rapprres[1]**2, (0.1/np.sqrt(3))**2])
    print(f"10*rapporto = {pr.Propagate(prodmean, prodcov, 2, prod, x)}")


    print("\nParte 1, Fabry")
    avg = 1/2*(x[0]+x[1])
    intextmeasures = [[3.5, 5.5, 6.8, 8], [4, 5.8, 7.3, 8.4]]
    avgcov = np.diag([(0.1/np.sqrt(3))**2]*2)
    avgarr = np.array([]).reshape(0,2)
    for i,j in zip(*intextmeasures):
        avgmean = [i, j]
        avgres = pr.Propagate(avgmean, avgcov, 2, avg, x)
        avgarr = np.vstack([avgarr, avgres])
    avgarr = np.vstack([avgarr, [9.2, 0.1/np.sqrt(3)], [10.1, 0.1/np.sqrt(3)]]).T
    costh = x[0]/pr.smp.sqrt(x[0]**2+x[1]**2/4)
    k = 10
    ang = pr.smp.acos(y[0])
    for i,j in zip(*avgarr):
        thmean = [207.4, i]
        thcov = np.diag([(0.1/np.sqrt(3))**2, j**2])
        thres = pr.Propagate(thmean, thcov, 2, costh, x)
        print(f"cos theta {k} = {thres}")
        angmean = [thres[0]]
        angcov = np.diag([thres[1]])
        print(f"theta {k} = {pr.Propagate(angmean, angcov, 1, ang, y)}")
        k -= 1

    print("\nParte 2, Fabry")
    thmean = [207.4, 4.6]
    thcov = np.diag([(0.1/np.sqrt(3))**2]*2)
    thres = pr.Propagate(thmean, thcov, 2, costh, x)
    print(f"cos theta = {thres}")

    deltad=x[0]*632.8/(2*x[1])
    dmean=[29.5,thres[0]]
    dcov=np.diag([0.453382350291181**2, thres[1]**2])
    deltad=x[0]*632.8e-9/(2*x[1])
    print(f"Delta d = {pr.Propagate(dmean,dcov,2,deltad,x)}")

    print("\nParte 3, Michelson")
    avgmean = [5.6, 4.4]
    avgcov = np.diag([(0.1/np.sqrt(3))**2]*2)
    avgres = pr.Propagate(avgmean, avgcov, 2, avg, x)

    costh = x[0]/pr.smp.sqrt(x[0]**2+x[1]**2/4)
    thmean = [196.7, avgres[0]]
    thcov = np.diag([(0.1/np.sqrt(3))**2, avgres[1]**2])
    thres = pr.Propagate(thmean, thcov, 2, costh, x)
    print(f"cos theta = {thres}")

    deltad=x[0]*632.8/(2*x[1])
    dmean=[29.4, thres[0]]
    dcov=np.diag([0.163299316185545**2, thres[1]**2])
    deltad=x[0]*632.8e-9/(2*x[1])
    print(f"Delta d = {pr.Propagate(dmean,dcov,2,deltad,x)}")

    print("\nParte 5, Vetro")
    costh = pr.smp.cos(y[0]*pr.smp.pi/180)
    thmean = [10]
    thcov = np.diag([(0.1/np.sqrt(3))**2])
    thres = pr.Propagate(thmean, thcov, 1, costh, y)
    print(f"cos theta = [{thres[0].evalf()}, {thres[1]}")

    d = 6e-3
    lam = 632.8e-9
    n=(2*d-x[0]*lam)*(1-x[1])/(2*d*(1-x[1])-x[0]*lam)
    nmean = [109.4, thres[0]]
    ncov = np.diag([1.0295630140987**2, thres[1]**2])
    nres = pr.Propagate(nmean, ncov, 2, n, x)
    print(f"n vetro = [{nres[0].evalf()}, {nres[1]}]")

    print("\nParte 7, righello")
    costh = x[0]/pr.smp.sqrt(x[0]**2+x[1]**2/4)
    thmean = [317, 46.5]
    thcov = np.diag([(0.5/np.sqrt(3))**2, (0.2/np.sqrt(3))**2])
    thres = pr.Propagate(thmean, thcov, 2, costh, x)
    print(f"cos theta incidente = {thres}")

    dists = [[49, 51.6, 53.8, 43.5, 40, 35.6], [(0.2/np.sqrt(3))**2]*6]
    distsoprhoriz = -x[0]/2.+x[1]
    disthoriz = np.array([]).reshape(0,2)
    for i, j in zip(*dists):
        distmean = [46.5, i]
        distcov = np.diag([(0.2/np.sqrt(3))**2, j**2])
        distres = pr.Propagate(distmean, distcov, 2, distsoprhoriz, x)
        disthoriz = np.vstack([disthoriz, distres])
    disthoriz = disthoriz.T

    k = a = 1
    costh = x[0]/pr.smp.sqrt(x[0]**2+x[1]**2)
    for i,j in zip(*disthoriz):
        thmean = [317, i]
        thcov = np.diag([(0.5/np.sqrt(3))**2, j**2])
        thres = pr.Propagate(thmean, thcov, 2, costh, x)
        print(f"cos theta {k} = {thres}")
        if k > 2:
            k = 0
            a = -1
        k += a

if __name__ == '__main__':
    main()