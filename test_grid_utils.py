from grid_utils import *
import pylab as plt

def test_meshgrid_pairs():
    phi = np.linspace(0, 2*np.pi, 361)
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    gp = meshgrid_pairs(theta, phi)

    assert gp.shape == (361, 181, 2)

def test_ungrid():
    phi = np.linspace(0, 2*np.pi, 361)
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    d = np.outer(np.sin(phi), np.sin(theta))
    xyz = ungrid(theta, phi, d)
    assert np.max(xyz[:, 0]) <= np.max(theta)  # X is theta
    assert np.max(xyz[:, 1]) <= np.max(phi)    # Y is phi
    assert xyz.shape == (361*181, 3)

def test_grid_xyz():
    phi = np.linspace(0, 2*np.pi, 361)
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    d = np.outer(np.sin(phi), np.sin(theta))
    xyz = ungrid(theta, phi, d)

    d_reconstructed = grid_xyz(xyz, 181, 361)

    plt.subplot(121)
    plt.imshow(d_reconstructed)
    plt.subplot(122)
    plt.imshow(d)
    plt.show()

    assert d_reconstructed.shape == d.shape

    assert np.allclose(d_reconstructed, d)

def test_scatter():
    a = np.array([1,2,3,4,5])
    idx = np.array([0,2,1,4,3])
    b = np.zeros_like(a)
    scatter(idx, a, b)
    assert np.allclose(b, [1,3,2,5,4])
    
    a = np.array([[1,2], [3,4]])
    idx = np.array([1,0])
    b = np.zeros_like(a)
    scatter(idx, a, b)
    assert np.allclose(b, np.array([[3,4], [1,2]]))
    
def test_gather():
    a = np.array([1,2,3,4,5])
    idx = np.array([0,2,1,4,3])
    b = gather(idx, a)
    assert np.allclose(b, [1,3,2,5,4])
    
    a = np.array([[1,2], [3,4]])
    idx = np.array([1,0])
    b = gather(idx, a)
    assert np.allclose(b, np.array([[3,4], [1,2]]))   

if __name__ == "__main__":
    test_scatter()
    test_gather()
    test_meshgrid_pairs()
    test_ungrid()
    test_grid_xyz()