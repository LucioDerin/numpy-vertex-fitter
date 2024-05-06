"""Script that produces the plots in ./images"""

from puma import Histogram, HistogramPlot
from puma.utils import get_good_colours
import numpy as np


def plot_Lxy_comparison(
    GN2_Lxy, perfect_tracksel_Lxy, SV1_Lxy, MCtruth_Lxy, title, filename
):
    """Plot the Lxy comparison histogram"""
    GN2_hist = Histogram(
        GN2_Lxy, label="Fit with GN2 track selection", histtype="step", alpha=1
    )
    ptracksel_hist = Histogram(
        perfect_tracksel_Lxy,
        label="Fit with perfect track selection",
        histtype="step",
        alpha=1,
    )
    SV1_hist = Histogram(SV1_Lxy, label="SV1", histtype="step", alpha=1)
    MCtruth_hist = Histogram(MCtruth_Lxy, label="MC truth", histtype="step", alpha=1)

    xrange = (0, 40)
    histogram_plotter = HistogramPlot(
        ylabel="Normalized Counts",
        xlabel=r"$L_{xy} [mm]$",
        logy=True,
        bins=xrange[1] - xrange[0],
        bins_range=xrange,
        norm=True,
        atlas_first_tag="Simulation Internal",
        atlas_second_tag=r"RUN3 $t\bar{t}$",
        figsize=(6, 5),
        n_ratio_panels=1,
    )
    histogram_plotter.title = title
    histogram_plotter.add(MCtruth_hist, reference=True)
    histogram_plotter.add(GN2_hist, reference=False)
    histogram_plotter.add(ptracksel_hist, reference=False)
    histogram_plotter.add(SV1_hist, reference=False)

    histogram_plotter.draw()
    filename = "images/" + filename
    histogram_plotter.savefig(filename, dpi=200, transparent=False)

def plot_residuals_comparison(
    GN2_Lxy, perfect_tracksel_Lxy, SV1_Lxy, MCtruth_Lxy, title, filename
):
    colors = get_good_colours()[1:]
    """Plot the Lxy comparison histogram"""
    GN2_hist = Histogram(
        (GN2_Lxy - MCtruth_Lxy), label="Fit with GN2 track selection", histtype="step", alpha=1,
        colour=colors[0],
    )
    ptracksel_hist = Histogram(
        (perfect_tracksel_Lxy-MCtruth_Lxy),
        label="Fit with perfect track selection",
        histtype="step",
        alpha=1,
        colour=colors[1]
    )
    SV1_hist = Histogram((SV1_Lxy - MCtruth_Lxy), label="SV1", histtype="step", alpha=1,colour=colors[2])

    xrange = (-40, 40)
    histogram_plotter = HistogramPlot(
        ylabel="Normalized Counts",
        xlabel=r"$fit_{L_{xy}} - MC_{L_{xy}} [mm]$",
        logy=True,
        bins=xrange[1] - xrange[0],
        bins_range=xrange,
        norm=True,
        atlas_first_tag="Simulation Internal",
        atlas_second_tag=r"RUN3 $t\bar{t}$",
        figsize=(6, 5),
        n_ratio_panels=1,
    )
    histogram_plotter.title = title
    histogram_plotter.add(GN2_hist, reference=False)
    histogram_plotter.add(ptracksel_hist, reference=True)
    histogram_plotter.add(SV1_hist, reference=False)

    histogram_plotter.draw()
    filename = "images/" + filename
    histogram_plotter.savefig(filename, dpi=200, transparent=False)


def plot_chi2_comparison(perfect_tracksel_chi2, GN2_chi2, title, filename):
    """Plot the chi2 comparison histogram"""
    GN2_hist = Histogram(
        GN2_chi2, label=r"$\chi^2$ with GN2 track selection", histtype="step", alpha=1
    )
    ptracksel_hist = Histogram(
        perfect_tracksel_chi2,
        label=r"$\chi^2$ with perfect track selection",
        histtype="step",
        alpha=1,
    )
    xrange = (0, 5)
    histogram_plotter = HistogramPlot(
        ylabel="Normalized Counts",
        xlabel=r"$\chi^2$",
        bins=(xrange[1] - xrange[0]) * 10,
        bins_range=xrange,
        logy=True,
        norm=True,
        atlas_first_tag="Simulation Internal",
        atlas_second_tag=r"RUN3 $t\bar{t}$",
        figsize=(6, 5),
        n_ratio_panels=1,
    )
    histogram_plotter.title = title
    histogram_plotter.add(ptracksel_hist, reference=True)
    histogram_plotter.add(GN2_hist, reference=False)
    histogram_plotter.draw()
    filename = "images/" + filename
    histogram_plotter.savefig(filename, dpi=200, transparent=False)


if __name__ == "__main__":
    GN2_Lxy = []
    GN2_chi2 = []
    perfect_tracksel_Lxy = []
    perfect_tracksel_chi2 = []
    SV1_Lxy = []
    MCtruth_Lxy = []
    jet_flav = []

    # Reading fit results
    with open("fit_results.dat", "r") as ifile:
        ifile.readline()
        for line in ifile:
            glxy, plxy, sv1lxy, tlxy, flav, TruthChi2, GN2chi2 = [
                float(i) for i in line.split(" ")
            ]
            GN2_Lxy.append(glxy)
            GN2_chi2.append(GN2chi2)
            perfect_tracksel_Lxy.append(plxy)
            perfect_tracksel_chi2.append(TruthChi2)
            SV1_Lxy.append(sv1lxy)
            MCtruth_Lxy.append(tlxy)
            jet_flav.append(flav)

    GN2_Lxy = np.array(GN2_Lxy)
    GN2_chi2 = np.array(GN2_chi2)
    perfect_tracksel_Lxy = np.array(perfect_tracksel_Lxy)
    perfect_tracksel_chi2 = np.array(perfect_tracksel_chi2)
    SV1_Lxy = np.array(SV1_Lxy)
    MCtruth_Lxy = np.array(MCtruth_Lxy)
    jet_flav = np.array(jet_flav)

    # Plots inclusive flavour
    plot_Lxy_comparison(
        GN2_Lxy,
        perfect_tracksel_Lxy,
        SV1_Lxy,
        MCtruth_Lxy,
        title=r"$L_{xy}$ comparison, $b/c$-jets",
        filename="LxyComparison_inclusive.jpg",
    )
    plot_residuals_comparison(
        GN2_Lxy,
        perfect_tracksel_Lxy,
        SV1_Lxy,
        MCtruth_Lxy,
        title=r"Residuals comparison, $b/c$-jets",
        filename="ResidualsComparison_inclusive.jpg",
    )
    plot_chi2_comparison(
        perfect_tracksel_chi2,
        GN2_chi2,
        title=r"$\chi^2$ comparison, $b/c$-jets",
        filename="chi2_inclusive.jpg",
    )

    # Plots only on c-jets
    filter_flav = 4
    mask = jet_flav == filter_flav
    GN2_Lxy_cjets = GN2_Lxy[mask]
    GN2_chi2_cjets = GN2_chi2[mask]
    perfect_tracksel_Lxy_cjets = perfect_tracksel_Lxy[mask]
    perfect_tracksel_chi2_cjets = perfect_tracksel_chi2[mask]
    SV1_Lxy_cjets = SV1_Lxy[mask]
    MCtruth_Lxy_cjets = MCtruth_Lxy[mask]
    plot_Lxy_comparison(
        GN2_Lxy_cjets,
        perfect_tracksel_Lxy_cjets,
        SV1_Lxy_cjets,
        MCtruth_Lxy_cjets,
        title=r"$L_{xy}$ comparison, $c$-jets only",
        filename="LxyComparison_cjets.jpg",
    )
    plot_residuals_comparison(
        GN2_Lxy_cjets,
        perfect_tracksel_Lxy_cjets,
        SV1_Lxy_cjets,
        MCtruth_Lxy_cjets,
        title=r"Residuals comparison, $c$-jets only",
        filename="ResidualsComparison_cjets.jpg",
    )
    plot_chi2_comparison(
        perfect_tracksel_chi2_cjets,
        GN2_chi2_cjets,
        title=r"$\chi^2$ comparison, $c$-jets only",
        filename="chi2_cjets.jpg",
    )

    # Plots only on b-jets
    filter_flav = 5
    mask = jet_flav == filter_flav
    GN2_Lxy_bjets = GN2_Lxy[mask]
    GN2_chi2_bjets = GN2_chi2[mask]
    perfect_tracksel_Lxy_bjets = perfect_tracksel_Lxy[mask]
    perfect_tracksel_chi2_bjets = perfect_tracksel_chi2[mask]
    SV1_Lxy_bjets = SV1_Lxy[mask]
    MCtruth_Lxy_bjets = MCtruth_Lxy[mask]
    plot_Lxy_comparison(
        GN2_Lxy_bjets,
        perfect_tracksel_Lxy_bjets,
        SV1_Lxy_bjets,
        MCtruth_Lxy_bjets,
        title=r"$L_{xy}$ comparison, $b$-jets only",
        filename="LxyComparison_bjets.jpg",
    )
    plot_residuals_comparison(
        GN2_Lxy_bjets,
        perfect_tracksel_Lxy_bjets,
        SV1_Lxy_bjets,
        MCtruth_Lxy_bjets,
        title=r"Residuals comparison, $b$-jets only",
        filename="ResidualsComparison_bjets.jpg",
    )
    plot_chi2_comparison(
        perfect_tracksel_chi2_bjets,
        GN2_chi2_bjets,
        title=r"$\chi^2$ comparison, $b$-jets only",
        filename="chi2_bjets.jpg",
    )
