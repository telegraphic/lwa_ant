from lwa_ant.lwa_observer import *
import time

def test_radec():

    elev = np.linspace(-np.pi/2, np.pi/2, 100000)
    az   = np.linspace(0, 2*np.pi, 100000)

    lwa_obs = LwaObserver()

    t1 = time.time()
    ra_decs = []
    for ii in range(len(elev)):
        a, e = az[ii], elev[ii]
        r, d = lwa_obs.radec_of(a, e)
        ra_decs.append((r, d))
    t2 = time.time()
    print "Time loop:       %2.2fs" % (t2 - t1)

    alt = elev
    t1 = time.time()
    ra_decs2 = lwa_obs.altaz2radec(alt, az)
    t2 = time.time()
    print "Time vectorized: %2.2fs" % (t2 - t1)

    for ii in range(len(ra_decs)):
        ra1, dec1 = float(ra_decs[ii][0]), float(ra_decs[ii][1])
        ra2, dec2 = ra_decs2[ii, 0], ra_decs2[ii, 1]
        assert ra1 == ra2
        assert dec1 == dec2

if __name__ == "__main__":
    test_radec()