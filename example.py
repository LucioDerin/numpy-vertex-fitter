# Modules import
from modules.ImportH5 import importH5
from modules.singleVertexFitter import singleVertexFitter_straightTracks as SVFs

# Python import
import numpy as np
import os
import sys

N = int(1e4)

filepath = None
# searching for the H5 database
for file in os.listdir():
    filename, ext = os.path.splitext(file)
    if ext == ".h5":
        filepath = file
        print(f"Found h5 database: {filepath}")

# If no h5 file, break
if filepath is None:
    sys.exit(1)

print("Importing jets.")
jets = importH5(
    filepath=filepath,
    Nevents=N,
    onlySV1=False,
    customProperties=["HadronConeExclTruthLabelLxy"],
)

# Single Secondary Vertex Fitter with straight line approximation
svfs = SVFs(eps=1e-6, maxIter=1e3)

# Number of successfully fitted jets
nFittedJets = 0
# Vertexes fitted with perfect track selection
perfect_tracksel_vertexes = []
# Chi2 of the perfect track selection fits
perfect_tracksel_chi2 = []
# Vertexes fitted with GN2 track selection
GN2_tracksel_vertexes = []
# Chi2 of the GN2 track selection fits
GN2_tracksel_chi2 = []
# Lxy by SV1
SV1_Lxy = []
# Montecarlo truth Lxy
MCtruth_Lxy = []
# Jet flavour label
jet_flavour = []

# Wether to fit or not light jets
filterLightJets = True

# Fitting jets
print("Begin fitting...")
for i, j in enumerate(jets):
    # Status
    if i % 50 == 0:
        print("Fitting jet", i, "out of", N, end="\r")

    # If enabled, skip light jets
    if j.properties["HadronConeExclTruthLabelID"] not in [4, 5] and filterLightJets:
        continue

    # Read track origins of the current jet and store the tracks in a list
    truth_track_origin = []
    GN2_track_origins = []
    tracks = []
    for t in j.allTracks:
        tracks.append(t)
        truth_track_origin.append(t.truthOriginLabel)
        GN2_track_origins.append(t.gn2Origin)

    truth_track_origin = np.array(truth_track_origin)
    GN2_track_origins = np.array(GN2_track_origins)
    tracks = np.array(tracks)

    # Track Selection with GN2
    maskGN2 = (
        (GN2_track_origins == "FromB")
        | (GN2_track_origins == "FromBC")
        | (GN2_track_origins == "FromC")
    )

    # Track Selection with MC truth origin
    maskTruth = (
        (truth_track_origin == "FromB")
        | (truth_track_origin == "FromBC")
        | (truth_track_origin == "FromC")
    )

    # tracks list masked with perfect track selection
    ptracksel_tracks = tracks[maskTruth]
    # tracks list masked with GN2 track selection
    GN2tracksel_tracks = tracks[maskGN2]

    # If no tracks are left, skip the jet
    if len(ptracksel_tracks) == 0 or len(GN2tracksel_tracks) == 0:
        continue

    # perfect tracksel fit
    vertexTruth, chi2Truth = svfs.fit(ptracksel_tracks.tolist())

    # GN2 tracksel fit
    vertexGN2, chi2GN2 = svfs.fit(GN2tracksel_tracks.tolist())

    # Saving results for this jet
    SV1_Lxy.append(j.properties["SV1_Lxy"])
    MCtruth_Lxy.append(j.properties["HadronConeExclTruthLabelLxy"])
    perfect_tracksel_vertexes.append(vertexTruth)
    perfect_tracksel_chi2.append(chi2Truth)
    GN2_tracksel_vertexes.append(vertexGN2)
    GN2_tracksel_chi2.append(chi2GN2)
    jet_flavour.append(j.properties["HadronConeExclTruthLabelID"])
    nFittedJets += 1

perfect_tracksel_vertexes = np.array(perfect_tracksel_vertexes)
GN2_tracksel_vertexes = np.array(GN2_tracksel_vertexes)
SV1_Lxy = np.array(SV1_Lxy)
MCtruth_Lxy = np.array(MCtruth_Lxy)
jet_flavour = np.array(jet_flavour)

# Calculating Lxy of the fitted vertex
# (for the coordinate system see H5Track docs)
GN2_Lxy = []
perfect_tracksel_Lxy = []
for i in range(nFittedJets):
    GN2_Lxy.append(
        np.linalg.norm([0, GN2_tracksel_vertexes[i, 1], GN2_tracksel_vertexes[i, 2]])
    )
    perfect_tracksel_Lxy.append(
        np.linalg.norm(
            [0, perfect_tracksel_vertexes[i, 1], perfect_tracksel_vertexes[i, 2]]
        )
    )

# Saving results
with open("fit_results.dat", "w") as ofile:
    print(
        "GN2_tracksel_Lxy perfect_tracksel_Lxy SV1_Lxy HadronConeExclTruthLabelLxy HadronConeExclTruthLabelID Truth_Chi2 GN2_chi2",
        file=ofile,
    )
    for i in range(nFittedJets):
        print(
            GN2_Lxy[i],
            perfect_tracksel_Lxy[i],
            SV1_Lxy[i],
            MCtruth_Lxy[i],
            jet_flavour[i],
            perfect_tracksel_chi2[i],
            GN2_tracksel_chi2[i],
            file=ofile,
        )
