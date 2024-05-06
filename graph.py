from puma import Histogram, HistogramPlot
from puma.utils import get_good_colours

import numpy as np

GN2_Lxy = []
perfect_tracksel_Lxy = []
SV1_Lxy = []
MCtruth_Lxy = []
jet_flav = []

with open("fit_results.dat", "r") as ifile:
    ifile.readline()
    for line in ifile:
        glxy, plxy, sv1lxy, tlxy, flav = [float(i) for i in line.split(" ")]
        GN2_Lxy.append(glxy)
        perfect_tracksel_Lxy.append(plxy)
        SV1_Lxy.append(sv1lxy)
        MCtruth_Lxy.append(tlxy)
        jet_flav.append(flav)

GN2_Lxy = np.array(GN2_Lxy)
perfect_tracksel_Lxy = np.array(perfect_tracksel_Lxy)
SV1_Lxy = np.array(SV1_Lxy)
MCtruth_Lxy = np.array(MCtruth_Lxy)
jet_flav = np.array(jet_flav)

filter_flav = 5
if filter_flav is not None:
    mask = (jet_flav==filter_flav)
    GN2_Lxy = GN2_Lxy[mask]
    perfect_tracksel_Lxy = perfect_tracksel_Lxy[mask]
    SV1_Lxy = SV1_Lxy[mask]
    MCtruth_hist = MCtruth_Lxy[mask]


GN2_hist = Histogram(GN2_Lxy,label="Fit with GN2 track selection", histtype="step", alpha=1)
ptracksel_hist = Histogram(perfect_tracksel_Lxy,label="Fit with perfect track selection", histtype="step", alpha=1)
SV1_hist = Histogram(SV1_Lxy,label="SV1", histtype="step", alpha=1)
MCtruth_hist = Histogram(MCtruth_Lxy,label="MC truth", histtype="step", alpha=1)

xrange = (0,40)
histogram_plotter = HistogramPlot(
    ylabel="Normalized Counts",
    xlabel=r"$L_{xy} [mm]$",
    logy=True,
    bins=xrange[1]-xrange[0],
    bins_range=xrange,
    norm=True,
    atlas_first_tag="Simulation Internal",
    atlas_second_tag=r"RUN3 $t\bar{t}$",
    figsize=(6, 5),
    n_ratio_panels=1,
)

histogram_plotter.add(MCtruth_hist, reference=True)
histogram_plotter.add(GN2_hist, reference=False)
histogram_plotter.add(ptracksel_hist, reference=False)
histogram_plotter.add(SV1_hist, reference=False)

histogram_plotter.draw()
filename = "images/LxyComparison.jpg" if filter_flav is None else f"images/LxyComparison_flav_{filter_flav}.jpg"
histogram_plotter.savefig(filename, dpi=200, transparent=False)