from lwa_ant.lwa_antenna import *

def test_generate():
    lwa = LwaBeamPattern()
    
    bp = lwa.generate(30)
    print bp.shape, lwa.phi.shape, lwa.theta.shape
    
    bp2 = lwa.generate(30, 181, 361)
    
    assert bp.shape == bp2.shape
    assert np.allclose(bp, bp2)
    
    plt.figure("Compare beam slices")
    plt.subplot(211)
    plt.plot(bp[0])
    plt.plot(bp2[0])
    plt.subplot(212)
    plt.plot(bp[90] - bp2[90])
    
    plt.show()

def test_view():
    lwa = LwaBeamPattern()
    lwa.generate(80)
    lwa.view()

def test_data():
    lwa = LwaBeamPattern()
    plt.imshow(lwa.data[0])
    plt.show()
    
    print lwa.data.shape

if __name__ == "__main__":
    test_generate()
    test_view()
    test_data()
    