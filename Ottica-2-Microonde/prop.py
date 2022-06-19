import sys
import numpy as np

sys.path.append('../Scripts/utils/')
import propagate as pr

def main():
    x = pr.GenerateXs(2)

    print("\nRifrazione:")
    #n1 = pr.smp.sin(x[0]*np.pi/180)/pr.smp.sin(x[1]*np.pi/180)
    #n1mean = [42.6666666666667, 30]
    #n1cov = np.diag([0.333333333333333**2, 2**2])
    #print(f"n1 = {pr.Propagate(n1mean, n1cov, 2, n1, x)}")
    n1 = pr.smp.sin((x[1]+180-x[0])*np.pi/180)/pr.smp.sin(x[1]*np.pi/180)
    n1mean = [167.333333333333, 22]
    n1cov = np.diag([0.333333333333333**2, 0**2])
    print(f"n1 = {pr.Propagate(n1mean, n1cov, 2, n1, x)}")

    print("\nBragg:")
    dist = 2.85*x[1]/(2*pr.smp.sin(x[0]*np.pi/180))
    distcov = np.diag([(1/np.sqrt(3))**2, 0])
    distmean = [22, 1]
    print(f"dist 22° = {pr.Propagate(distmean, distcov, 2, dist, x)}")
    distmean = [47, 2]
    print(f"dist 47° = {pr.Propagate(distmean, distcov, 2, dist, x)}")

    print("\nDouble slit:")
    wave = 7.6*pr.smp.sin(x[0]*np.pi/180)/x[1]
    wavecov = np.diag([(1/np.sqrt(3))**2, 0])
    wavemean = [22, 1]
    print(f"wave 22° = {pr.Propagate(wavemean, wavecov, 2, wave, x)}")
    wavemean = [48, 2]
    print(f"wave 48° = {pr.Propagate(wavemean, wavecov, 2, wave, x)}")

    print("\nLloyd:")
    x = pr.GenerateXs(6)
    #dist_centr_58 = 5.96
    #postrasm = 12
    #h1mean = 78.06
    #dst_centro_specchio = 7.256
    #posricv = 112
    #dist_centr_69 = 7.25

    #poscentro_58 = dist_centr_58+58
    #poscentro_69 = 69-dist_centr_69
    #posricv_real = posricv+poscentro_58-poscentro_69
    #diff = h1mean-69+dst_centro_specchio
    #d1 = poscentro_58-postrasm
    #d2 = posricv_real-poscentro_58
    #BC = sqrt(d2**2+diff**2)
    #BC = sqrt((posricv-(69-dist_centr_69))**2+(h1mean-69+dst_centro_specchio)**2)
    #AB = sqrt(d1**2+diff**2)
    #AB = sqrt((dist_centr_58+58-postrasm)**2+(h1mean-69+dst_centro_specchio)**2)
    #leng = AB+BC

    # x0 = distanza centro goniometro e segno 58cm sull'asta del trasmettitore
    # x1 = posizione trasmettitore
    # x2 = distanza specchio media letta sul metro (non effettiva)
    # x3 = distanza centro goniometro e segno 69cm sull'asta dello specchio
    # x4 = posizione ricevitore
    # x5 = distanza centro gonomiometro e segno 69cm sull'asta del ricevitore
    length = pr.smp.sqrt((x[0]+58-x[1])**2+(x[2]-69+x[3])**2)+pr.smp.sqrt((x[4]-(69-x[5]))**2+(x[2]-69+x[3])**2)
    mean = [5.96, 12, 78.06, 7.256, 112, 7.25]
    cov = np.diag([(0.01/np.sqrt(3))**2, (0.1/np.sqrt(3))**2, 0.102956301409869**2, (0.002/np.sqrt(3))**2, (0.1/np.sqrt(3))**2, (0.01/np.sqrt(3))**2])
    res1 = pr.Propagate(mean, cov, 6, length, x)
    print(f"AB+BC 1: {res1}")

    mean[2] = 82.14
    cov[2, 2] = 0.116619037896907**2
    res2 = pr.Propagate(mean, cov, 6, length, x)
    print(f"AB+BC 2: {res2}")

    x = pr.GenerateXs(2)
    wavelength = x[1]-x[0]
    wmean = [res1[0], res2[0]]
    wcov = np.diag([res1[1]**2, res2[1]**2])
    print(f"wavelength: {pr.Propagate(wmean, wcov, 2, wavelength, x)}")

if __name__ == '__main__':
    main()