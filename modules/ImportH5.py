# Modules import
from modules.containers import H5Track
from modules.containers import JetContainer

# Python import
import h5py
import sys
import math as m
import numpy as np

def importH5(
    filepath="",
    Nevents: int = -1,
    customProperties: list = [],
    onlySV1: bool = True,
    straightTracks: bool = True,
):
    """Function that imports jets from an H5 file and store them in a list of JetContainer.

    Parameters
    ----------
    filepath : str, optional
        Path to the H5 file, by default ""
    Nevents : int, optional
        Number of jet to be read (does NOT correspond to the number of 
        imported jets if onlySV1==True), by default -1 which means import all dataset
    customProperties : list, optional
        other jet properties to import aside from the default of JetContainer;
        specify them by their key name in the H5 file, by default []
    onlySV1 : bool, optional
        Wether to filter only jets that have been fitted by SV1, by default True
    straightTracks : bool, optional
        Wether to use straight (True) or curved (False) tracks, by default True

    Returns
    -------
    list
        a list of JetContainer.
    """
    if onlySV1:
        print(
            "Warning: filtering out jets that have no SV1 fit. This might lead to fewer imported jets than specified in Nevents."
        )

    # File reading
    try:
        with h5py.File(filepath, "r") as h5Database:
            if Nevents != -1:
                jets = h5Database["jets"][:Nevents]
                rawTracks = (
                    h5Database["tracks_loose"][:Nevents]
                    if "tracks_loose" in h5Database.keys()
                    else h5Database["tracks"][:Nevents]
                )
            else:
                jets = h5Database["jets"][:]
                rawTracks = (
                    h5Database["tracks_loose"][:]
                    if "tracks_loose" in h5Database.keys()
                    else h5Database["tracks"][:]
                )

    except Exception as e:
        print("Error reading the file!")
        print(e)
        print("Ended exception")
        sys.exit(1)

    # Loop to import as many jets as possible (but < Nevents)
    importedJets = []
    # iterating over jets in the H5 file
    for jet, jetRawTracks in zip(jets, rawTracks):
        # Saving jet's variables based on fields' names
        jetKeys = [str(name) for name in jet.dtype.names]
        nTracks = (
            jet[jetKeys.index("n_tracks_loose")]
            if "n_tracks_loose" in jetKeys
            else jet[jetKeys.index("n_tracks")]
        )
        SV1_L3d = jet[jetKeys.index("SV1_L3d")]

        # filtering only jets fitted by SV1, if onlySV1==True
        if onlySV1 and m.isnan(SV1_L3d):
            continue

        # Reading other jet's properties
        SV1_Lxy = jet[jetKeys.index("SV1_Lxy")]
        jetEta = jet[jetKeys.index("eta")]
        jetPhi = jet[jetKeys.index("phi")] if "phi" in jetKeys else 0
        primaryVertexDetectorZ = (
            jet[jetKeys.index("primaryVertexDetectorZ")]
            if "primaryVertexDetectorZ" in jetKeys
            else 0
        )
        jetPt = jet[jetKeys.index("pt")]
        HadronConeExclTruthLabelID = jet[jetKeys.index("HadronConeExclTruthLabelID")]

        # if specified, import custom properties of the jet:
        customPropertiesDict = {}
        if len(customProperties) != 0:
            for p in customProperties:
                customPropertiesDict[p] = jet[jetKeys.index(p)]

        tracksSV1 = []
        tracksNoSV1 = []
        # Iterating over non zero-padded tracks of the jet
        for rawTrack in jetRawTracks[:nTracks]:
            # If using straight tracks
            if straightTracks:
                # Saving track's fields' name
                trackKeys = [str(name) for name in rawTrack.dtype.names]
                pt = rawTrack[trackKeys.index("pt")]
                # Reading track's variables
                eta = rawTrack[trackKeys.index("eta")]
                # Rotating the track so that the jet is displayed vertically
                dphi = rawTrack[trackKeys.index("dphi")] + m.pi / 2.0
                IP3D_signed_d0 = rawTrack[trackKeys.index("IP3D_signed_d0")]
                z0RelativeToBeamspot = rawTrack[trackKeys.index("z0RelativeToBeamspot")]
                truthOriginLabel = rawTrack[trackKeys.index("ftagTruthOriginLabel")]
                SV1VertexIndex = rawTrack[trackKeys.index("SV1VertexIndex")]

                # Track's error
                sigmaz0 = rawTrack[trackKeys.index("z0RelativeToBeamspotUncertainty")]
                sigmaPhi = rawTrack[trackKeys.index("phiUncertainty")]
                sigmaTheta = rawTrack[trackKeys.index("thetaUncertainty")]
                sigmaD0 = rawTrack[trackKeys.index("d0Uncertainty")]
                error = [sigmaTheta, sigmaPhi, sigmaD0, sigmaz0]

                # Track origin prediction
                originProbabilities = []
                originProbabilities.append(rawTrack[trackKeys.index("Pileup")])
                originProbabilities.append(rawTrack[trackKeys.index("Fake")])
                originProbabilities.append(rawTrack[trackKeys.index("Primary")])
                originProbabilities.append(rawTrack[trackKeys.index("FromB")])
                originProbabilities.append(rawTrack[trackKeys.index("FromBC")])
                originProbabilities.append(rawTrack[trackKeys.index("FromC")])
                originProbabilities.append(rawTrack[trackKeys.index("FromTau")])
                originProbabilities.append(rawTrack[trackKeys.index("OtherSecondary")])

                predictedOrigin = np.argmax(np.array(originProbabilities))

                # Creating the H5Track object
                track = H5Track(
                    pt,
                    eta,
                    dphi,
                    IP3D_signed_d0,
                    z0RelativeToBeamspot,
                    error,
                    truthOriginLabel,
                    predictedOrigin,
                    SV1VertexIndex,
                )

                # Saving the track in its respective list (selected by SV1 or not)
                if SV1VertexIndex == 0:
                    tracksSV1.append(track)
                else:
                    tracksNoSV1.append(track)
            else:
                # Implement charged tracks import
                # for non linear vertex fit
                pass

        # Filtering events that have a SV1 tracks list
        if onlySV1 and len(tracksSV1) == 0:
            continue

        # If jet is valid, create its jet container
        importedJet = JetContainer(
            tracksNoSV1,
            tracksSV1,
            nTracks=nTracks,
            SV1_L3d=SV1_L3d,
            SV1_Lxy=SV1_Lxy,
            eta=jetEta,
            phi=jetPhi,
            pt=jetPt,
            primaryVertexDetectorZ=primaryVertexDetectorZ,
            HadronConeExclTruthLabelID=HadronConeExclTruthLabelID,
            **customPropertiesDict,
        )
        # Store jet's container
        importedJets.append(importedJet)

    print("Successfully imported", len(importedJets), "jets.")

    # Return the list of JetContainer
    return importedJets
