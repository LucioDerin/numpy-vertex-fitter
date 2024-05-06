# Python imports
import numpy as np
import math as m

class JetContainer:
    """Class that contains the important elements of a jet in a vertex fitting perspective. 
    @self.tracksNoSV1: list of H5Tracks that are not selected by SV1;
    @self.tracksSV1: list of H5Tracks that are selected by SV1;
    @self.properties: properties of the jet as a whole (extracted from the h5); at the moment they are:
            - nTracks;
            - SV1_L3d;
            - SV1_Lxy;
            - eta;
            - phi;
            - pt;
            - primaryVertexDetectorZ;
            - HadronConeExclTruthLabelID;
        Keeping the naming coherent with H5 files.
    """
    def __init__(self,tracksNoSV1,tracksSV1,**jetProperties) -> None:
        self.properties = jetProperties
        self.tracksNoSV1 = tracksNoSV1
        self.tracksSV1 = tracksSV1
        self.allTracks = self.tracksNoSV1 + self.tracksSV1

class H5Track:
    '''
    Class that implements the representation of a track as a straight line with an origin and a versor.
    Note: coordinates are stored with the axes' names choice coherent to matplotlib 3D axes: right-handed
    with z as the vertical axis, which differs from LHC choice:
        LHC  -> MATPLOTLIB
        x    ->    y
        y    ->    z
        z    ->    x

    Public Members:
    @self.origin: np.array of shape (3,), (x,y,z) coordinates of the track's origin (perigee wrt beamline);
    @self.versor: np.array of shape (3,), (x,y,z) components of the track's versor
        with LHC choice of reference system they are (cos(theta),cos(phi)sin(theta),sin(phi)sin(theta));
    @self.truthOriginLabel: str indicating the track's provenience's MC truth;
    @self.SV1VertexIndex: str indicating the SV1 selection for the track's vertex ('From PV' or 'From SV');

    Constructors:
    __init__(eta, dphi, IP3D_signed_d0, z0RelativeTobeamspot): builds the track from its parameter in the H5 file;

    Public Methods:
    @evaluate(t): returns np.array of shape (3,), evaluates the parametric representation of the line 
        where t is the value of the line's parameter;

    Private Methods:
    @__theta(eta): converts the pseudorapidity into polar angle theta;
    '''

    # Dictionary to map the numeric value of the track's truthOriginDict to its meaning
    truthOriginDict = {
        -1:'ND',
        0: 'PU',
        1: 'Fake',
        2: 'Primary',
        3: 'FromB',
        4: 'FromBC',
        5: 'FromC',
        6: 'FromTau',
        7: 'OtherSecondary',
        8: 'ND'
    }

    # Private methods
    # eta -> theta conversion
    def __theta(self, eta):
        return 2*m.atan(m.exp(-eta))

    # Constructor
    def __init__(self, pt, eta, dphi, IP3D_signed_d0, z0RelativeToBeamspot, errors, truthOriginLabel=8, gn2Origin = 8, SV1VertexIndex=-2):
        '''
        Constructor of the class from h5 track's parameter.
        Parameters:
        @pt: track's pt;
        @eta: tracks's eta field in h5 file;
        @dphi: tracks's dphi field in h5 file;
        @IP3D_signed_d0: tracks's IP3D_signed_d0 field in h5 file;
        @z0RelativeToBeamspot: tracks's z0RelativeToBeamspot field in h5 file;
        @truthOriginLabel: tracks's truthOriginLabel field in h5 file;
        @SV1VertexIndex: tracks' SV1VertexIndex field in h5 file;
        '''
        self.pt = pt
        # eta -> theta conversion
        thetaT = self.__theta(eta)
        # versor
        xCern = m.sin(thetaT)*m.cos(dphi)
        yCern = m.sin(thetaT)*m.sin(dphi)
        zCern = m.cos(thetaT)
        self.versor = np.array([zCern, xCern, yCern])

        # origin
        phiP = dphi - m.pi/2 if IP3D_signed_d0 > 0 else dphi + m.pi/2
        y0 = abs(IP3D_signed_d0)*m.sin(phiP)
        x0 = abs(IP3D_signed_d0)*m.cos(phiP)
        self.origin = np.array([z0RelativeToBeamspot, x0, y0])
        self.IP3D_signed_d0 = IP3D_signed_d0
        # Errors
        sigmaTheta = errors[0]
        sigmaPhi = errors[1]
        sigmaD0 = errors[2]
        sigmaZ0 = errors[3]

        # versor error
        eVersor = np.array([((m.cos(thetaT)*m.cos(dphi)*sigmaTheta)**2) + ((m.sin(dphi)*m.sin(thetaT)*sigmaPhi)**2),
                                 ((m.cos(thetaT)*m.sin(dphi)*sigmaTheta)**2) + ((m.sin(thetaT)*m.cos(dphi)*sigmaPhi)**2),
                                 (m.sin(thetaT)*sigmaTheta)**2])
        # origin error
        eOrigin = np.array([sigmaZ0**2,
                            (m.cos(phiP)*sigmaD0)**2 + (m.sin(phiP)*IP3D_signed_d0*sigmaPhi)**2,
                            (m.sin(phiP)*sigmaD0)**2 + (m.cos(phiP)*IP3D_signed_d0*sigmaPhi)**2])
        
        eVersor = np.sqrt(eVersor)
        eOrigin = np.sqrt(eOrigin)

        # cov matrix
        # for now, I only have the diagonal
        self.covMat = np.diag(np.concatenate((eOrigin,eVersor)))

        # MC truth and GN2 label
        self.truthOriginLabel = self.truthOriginDict[truthOriginLabel]
        self.gn2Origin = self.truthOriginDict[gn2Origin]
        self.SV1VertexIndex = "From SV" if SV1VertexIndex==0 else "From PV"

    def evaluate(self, t):
        '''
        Evaluates the parametric representation of the track at the parameter value t.
        Parameters:
        @t: value of the parameter of the parametric representation of the line;
        Returns:
        @P: np.array of shape (3,), the track's point at parameter = t;
        '''
        P = self.origin + self.versor*t
        return P
    
    def print(self, ofile=None):
        if ofile != None:
            print(self.origin[0],self.origin[1],self.origin[2],self.versor[0],self.versor[1],self.versor[2],file=ofile)
            print(self.eOrigin[0],self.eOrigin[1],self.eOrigin[2],self.eVersor[0],self.eVersor[1],self.eVersor[2],file=ofile)
        else:
            print(self.origin[0],self.origin[1],self.origin[2],self.versor[0],self.versor[1],self.versor[2])
            print(self.eOrigin[0],self.eOrigin[1],self.eOrigin[2],self.eVersor[0],self.eVersor[1],self.eVersor[2])
        return