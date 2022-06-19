import sys
import numpy as np

sys.path.append('../Scripts/utils/')
import propagate as pr

def main():
    x = pr.GenerateXs(2)

    print("\nIndici di rifrazione prisma mercurio:")
    angprisma = 60
    torad = np.pi/360
    angs = [272.3333333333333333, 271.6666666666666666, 271.3333333333333333, 269.3333333333333333, 269]
    n = pr.smp.sin((x[0]-x[1]+angprisma)*torad)/pr.smp.sin(angprisma*torad)
    ncov = np.diag([1/(3*np.sqrt(3))**2]*2)
    i = 1
    nmean = [319.3333333333333333, 0]
    for j in angs:
        nmean[1] = j
        print(f"n {i} = {pr.Propagate(nmean, ncov, 2, n, x)}")
        i += 1

    print("\nIndici di rifrazione prisma ignoto:")
    angs = [273, 272.6666666666666666, 272, 271.3333333333333333, 271, 270.6666666666666666, 270]
    i = 1
    res = []
    nmean = [319.3333333333333333, 0]
    for j in angs:
        nmean[1] = j
        res.append(pr.Propagate(nmean, ncov, 2, n, x))
        print(f"n {i} = {res[i-1]}")
        i += 1

    print("\nLunghezze d'onda tubo spettrale:")
    x = pr.GenerateXs(3)
    lam = pr.smp.sqrt(x[1]/(x[2]-x[0]))
    lammean = [1.58933434040652, 0.881417892557338, 0]
    lamcov = np.array([[0.00347669143634916**2, -0.000264530216674947, 0],\
        [-0.000264530216674947, 0.0814341306040933**2, 0],\
        [0, 0, 0]])
    for j in range(len(angs)):
        lammean[2] = res[j][0]
        lamcov[2,2] = res[j][1]**2
        print(f"Lambda {j+1} = {pr.Propagate(lammean, lamcov, 3, lam, x)}")
    
    print("\nPasso reticolo 300:")
    lamb = 589.2937
    angs = [310, 299.45, 288.16666666666666, 275, 259.33333333333333333]
    d = x[0]*lamb/pr.smp.sin((x[1]-x[2])*np.pi/180)
    dmean = [0, 320, 0]
    dcov = np.diag([0, 1/(3*np.sqrt(3))**2, 1/(3*np.sqrt(3))**2])
    print("Destra:")
    for j in range(len(angs)):
        dmean[0] = j+1
        dmean[2] = angs[j]
        print(f"Passo {j+1} = {pr.Propagate(dmean, dcov, 3, d, x)}")
    angs = [329.666666666666666, 340.66666666666666666, 352, 364.66666666666666666]
    d = x[0]*lamb/pr.smp.sin((x[2]-x[1])*np.pi/180)
    print("Sinistra:")
    for j in range(len(angs)):
        dmean[0] = j+1
        dmean[2] = angs[j]
        print(f"Passo {j+1} = {pr.Propagate(dmean, dcov, 3, d, x)}")

    print("\nPasso reticolo 600:")
    lamb = 589.2937
    angs = [300.333333333333333, 277]
    d = x[0]*lamb/pr.smp.sin((x[1]-x[2])*np.pi/180)
    dmean = [0, 319.3333333333333333, 0]
    dcov = np.diag([0, 1/(3*np.sqrt(3))**2, 1/(3*np.sqrt(3))**2])
    print("Destra:")
    for j in range(len(angs)):
        dmean[0] = j+1
        dmean[2] = angs[j]
        print(f"Passo {j+1} = {pr.Propagate(dmean, dcov, 3, d, x)}")
    angs = [340.666666666666666, 369.33333333333333333]
    d = x[0]*lamb/pr.smp.sin((x[2]-x[1])*np.pi/180)
    print("Sinistra:")
    for j in range(len(angs)):
        dmean[0] = j+1
        dmean[2] = angs[j]
        print(f"Passo {j+1} = {pr.Propagate(dmean, dcov, 3, d, x)}")

    print("\nLunghezze d'onda tubo spettrale:")
    d = 1e-3/300
    angs = [311.666666666666, 311.33333333333333, 311, 310.33333333333333]
    lam = d*pr.smp.sin((x[0]-x[1])*np.pi/180)/x[2]
    lammean = [320.33333333333333333, 0, 1]
    lamcov = np.diag([1/(3*np.sqrt(3))**2, 1/(3*np.sqrt(3))**2, 0])
    print("Destra:")
    print("n = 1:")
    lammean[2] = 1
    for j in range(len(angs)):
        lammean[1] = angs[j]
        print(f"Lambda {j+1} = {pr.Propagate(lammean, lamcov, 3, lam, x)}")
    angs = [304, lammean[0], 301.66666666666666, 299.3333333333333333]
    print("n = 2:")
    lammean[2] = 2
    for j in range(len(angs)):
        lammean[1] = angs[j]
        print(f"Lambda {j+1} = {pr.Propagate(lammean, lamcov, 3, lam, x)}")
    print("Sinistra:")
    lam = d*pr.smp.sin((x[1]-x[0])*np.pi/180)/x[2]
    angs = [327.33333333333333, lammean[0], 328, 329.66666666666]
    print("n = 1:")
    lammean[2] = 1
    for j in range(len(angs)):
        lammean[1] = angs[j]
        print(f"Lambda {j+1} = {pr.Propagate(lammean, lamcov, 3, lam, x)}")
    angs = [335, 335.6666666666, 337, 340]
    print("n = 2:")
    lammean[2] = 2
    for j in range(len(angs)):
        lammean[1] = angs[j]
        print(f"Lambda {j+1} = {pr.Propagate(lammean, lamcov, 3, lam, x)}")

if __name__ == '__main__':
    main()