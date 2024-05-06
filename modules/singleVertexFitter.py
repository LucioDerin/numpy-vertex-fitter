import numpy as np


class singleVertexFitter_straightTracks:
    """Class that implements a Single Vertex Fitter on straight tracks based on least squares minimization."""

    def __init__(self, eps: float = 1e-8, maxIter: float = 1e2) -> None:
        """Constructor of the fitter.

        Parameters
        ----------
        eps : float, optional
            stability tolerance that, when reached, stops the fit, by default 1e-8
        maxIter : float, optional
            maximum number of minimization iteration, by default 1e2
        """
        self.eps = eps
        self.maxIter = maxIter

    def __Di(self, r, v, a):
        Di = np.cross((r - v), a) ** 2
        return Di

    def fit(self, tracks: list):
        """Functions that fit a single vertex on the tracks (H5Tracks)
        contained in the given list, assumed to be straight tracks.

        Parameters
        ----------
        tracks : list
            tracks (list of H5Tracks): list of tracks to be fitted

        Returns
        -------
        np.ndarray of shape (3,)
            coordinates of the fitted vertex [z,x,y].
        """
        # Handy variables
        Ntracks = len(tracks)
        iter = 0
        dv = 100
        v = np.zeros(3)

        # Initialize the vertex in the average origin of the tracks
        for t in tracks:
            v = v + (t.origin)
        v = v / (Ntracks)
        # Number of times the vertex has not been significantly moved
        nStuck = 0

        # Minimization loop
        while iter < self.maxIter:
            # Current jacobian
            J = []
            # Current derivative wrt tracks origin
            Dr = []
            # If no significant increment wrt the previous iteration
            if dv < self.eps:
                # increment the counter
                nStuck += 1
                # if maximum tries are reached, stop
                if nStuck > 5:
                    break
            # Else, resetting the counter
            else:
                nStuck = 0

            # Iterate over tracks to build the Jacobian of the least squares
            for t in tracks:
                # auxiliary variables defined in the literature
                di = t.origin - v
                eta = np.zeros((3, 3))
                ai = t.versor

                # Building eta
                for j in range(3):
                    for k in range(3):
                        eta[j, k] = ai[j] * di[k] - ai[k] * di[j]

                # Derivative wrt track's origin
                Dir = [
                    ai[1] * eta[1, 0] - ai[2] * eta[0, 2],
                    ai[2] * eta[2, 1] - ai[0] * eta[1, 0],
                    ai[0] * eta[0, 2] - ai[1] * eta[2, 1],
                ]

                # Derivative wrt track's versor
                Dia = [
                    di[2] * eta[0, 2] - di[1] * eta[1, 0],
                    di[0] * eta[1, 0] - di[2] * eta[2, 1],
                    di[1] * eta[2, 1] - di[0] * eta[0, 2],
                ]

                Dir = 2 / Ntracks * np.array(Dir)
                Dia = 2 / Ntracks * np.array(Dia)

                # Storing derivative
                Dr.append(Dir)
                # Building the Jacobian
                Ji = np.concatenate((Dir, Dia), axis=0)
                J.append(Ji)

            # Calculating tracks' weights
            sigma = []
            for i, Ji in enumerate(J):
                sigmai = Ji.T @ tracks[i].covMat @ Ji
                # stability for the inversion
                sigmai += 1e-9
                sigma.append(sigmai)

            # Calculating the gradient of the least squares
            gradS = np.zeros(3)
            for Dir, sigmai in zip(Dr, sigma):
                gradS = gradS - (sigmai**-1) * Dir

            # Calculating the Hessian of the least squares
            laplS = np.zeros((3, 3))
            # Iterating over each track to build the respective hessian
            for i, t in enumerate(tracks):
                ai = t.versor
                Hi = np.diag(
                    [
                        ai[1] ** 2 + ai[2] ** 2,
                        ai[0] ** 2 + ai[2] ** 2,
                        ai[0] ** 2 + ai[1] ** 2,
                    ]
                )
                Hi[0, 1] = Hi[1, 0] = -ai[0] * ai[1]
                Hi[0, 2] = Hi[2, 0] = -ai[0] * ai[2]
                Hi[1, 2] = Hi[2, 1] = -ai[1] * ai[2]
                laplS = laplS + 2 / Ntracks * (sigma[i] ** -1) * Hi

            # Storing old vertex
            vold = np.copy(v)
            # Updating the vertex
            v -= np.linalg.pinv(laplS) @ gradS
            # Calculating the increment
            dv = np.linalg.norm((vold - v))
            # Incrementing iteration counter
            iter += 1

        D = []
        for t in tracks:
            Di = self.__Di(t.origin, v, t.versor)
            D.append(np.sqrt(np.sum(Di**2)))
        chi2 = np.sum(np.array(D))

        # Returning vertex and chi2
        return np.array(v), chi2
