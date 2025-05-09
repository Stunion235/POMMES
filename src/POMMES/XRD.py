import numpy as np
import math
import readCIF
def interpolate(disordered,ordered,order):
    """
    Interpolate any property between the disordered and ordered structures.

    Parameters
    ----------
    disordered : array_like
        Value of the property in the disordered cell.
    ordered : array_like
        Value of the property in the ordered cell.
    order : float
        Fraction of order from 0 to 1.

    Returns
    -------
    array_like
        Interpolation of property at the given amount of order.

    """
    if order > 1 or order < 0:
        raise ValueError("Invalid order value")
    return order*ordered+(1-order)*disordered
def plane_spacing(hkl,abc):
    """
    Find the Plane spacing in Å of a set of planes.

    Parameters
    ----------
    hkl : list[int]
        Miller indices [h,k,l] of the diffraction plane.
    abc : list[float]
        Lattice constants [a,b,c].

    Returns
    -------
    float
        Plane spacing in Å for the family of planes.

    """
    return 1/(math.sqrt((hkl[0]/abc[0])**2 + (hkl[1]/abc[1])**2 + (hkl[2]/abc[2])**2))
def f_from_params(f,d):
    """
    Find the form factor from a set of form factor parameters.

    Parameters
    ----------
    f : list[float]
        List of form factor parameters a1,b1,a2,b2...c.
    d : float
        Plane spacing in Å.

    Returns
    -------
    float
        The form factor.

    """
    #Equivalent to ∑(a*e^(-b*(Q/(4π))^2))+c
    Q = (2*math.pi)/d
    while (np.shape(f)[0]%2)!=0:
        f=np.append(f,0)
    f=f.reshape((int(np.shape(f)[0]/2),2))
    return np.sum(f[:,0] * np.exp(-1 * f[:,1] * (Q/(4*math.pi))**2))
def form_factors(atoms,d):
    """
    Get the form factors for a list of atoms

    Parameters
    ----------
    atoms : list[str]
        List of atoms.
    d : float
        Plane spacing in Å.

    Returns
    -------
    list[float]
        List of form factors in the order given by `atoms`

    """
    factors=[]
    formcsv=open('../../data/atomic_formfactors.csv','r')
    forms=formcsv.readlines()
    for t in atoms:
        rownum=-1
        try:
            rownum=readCIF.line_index_of(forms, t)
        except IndexError:
            try:
                rownum=readCIF.line_index_of(forms, t[:-2])
            except IndexError:
                readCIF.warn_user(t+' is not a recognized atom/ion.', 'Continue, assuming form factor of 0?')
        if (rownum<1):
            factors.append(0)
        else:
            row=list(map(float,forms[rownum].split(',')[1:]))
            factors.append(f_from_params(row,d))
    formcsv.close()
    return(factors)
def structure_factor(formfactors,disordered_matrix,ordered_matrix,hkl,order):
    """
    Get the structure factor for a unit cell.

    Parameters
    ----------
    formfactors : list[float]
        List of form factors in the order given by the matrices.
    disordered : np.ndarray
        Matrix where columns 0, 1, and 2 are the x, y, and z fractional
        coordinates, and column 3 is the occupancy, of each unique site/atom
        combination in the disordered CIF in the order given by `sites`.
    ordered : np.ndarray
        Same as `disordered` but for the ordered CIF.
    hkl : list[int]
        Miller indices [h,k,l] of the diffraction plane.
    order : float
        Amount of order in the structure.

    Returns
    -------
    float
        The unit cell's structure factor.

    """
    factor=-1j*2*math.pi
    Fs=np.zeros_like(formfactors,dtype="complex")
    for i in range(len(formfactors)):
        Fs[i]=interpolate(disordered_matrix[i,-1],ordered_matrix[i,-1],order)*formfactors[i]*math.e**(factor*np.dot(hkl,interpolate(disordered_matrix[i,0:3],ordered_matrix[i,0:3],order)))
    return sum(Fs)
def mult(hkl):
    """
    Get the multiplicity of a peak (number of equivalent peaks)

    Parameters
    ----------
    hkl : list[int]
        Miller indices [h,k,l] of the diffraction plane.

    Returns
    -------
    int
        Multiplicity.

    """
    return 2**(len(list(filter(lambda a:a!=0,hkl))))
def intensity(S,d,wavelength):
    """
    Get the intensity of a peak given the structure factor and Plane spacing in Å.

    Parameters
    ----------
    S : float
        Unit cell's structure factor.
    d : float
        Plane spacing in Å.
    wavelength : float
        X-ray wavelength in Å.

    Returns
    -------
    float
        Expected intensity of the peak.

    """
    theta_rad = math.asin(wavelength/(2*d))
    twotheta_rad = theta_rad*2
    theta_M_rad = 22.65*(math.pi/180)
    return ((S**2) * (1+(math.cos(2*theta_M_rad)**4)*(math.cos(twotheta_rad)**2))/(math.sin(theta_rad)*math.cos(theta_rad)))