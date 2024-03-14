import numpy as np
import plotly.express as px
from scipy.interpolate import NearestNDInterpolator
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas

def plot_dose_map_contours(dose_map_to_plot,levels=3,**kwargs):

    dose_map_to_plot_sorted = dose_map_to_plot.sort_values(by=["longitudeTranslated","latitude"])

    (contour_longs, contour_lats) = np.meshgrid(dose_map_to_plot_sorted["longitudeTranslated"].unique(),
            dose_map_to_plot_sorted["latitude"].unique())

    interp = NearestNDInterpolator(list(zip(dose_map_to_plot_sorted["longitudeTranslated"], dose_map_to_plot_sorted["latitude"])),
                               dose_map_to_plot_sorted["edose"])

    contours = plt.contour(contour_longs,contour_lats,interp(contour_longs, contour_lats),
            levels=levels,linestyles="dashed",colors="black",zorder=1000,**kwargs)
    plt.clabel(contours, inline=True) #,fmt={"fontweight":"bold"})

def create_single_dose_map_plot_plt(heatmap_DF_to_Plot,
                                    hue_range = None, #(0,100), 
                                    heatmap_s = 63,
                                    edgecolor=None,
                                    dose_type = "edose",
                                    legend_label=r"Effective dose ($\mu Sv / hr$)",
                                    palette="Spectral_r",
                                    plot_longitude_east=False,
                                    plot_colorbar=True):

    if not (heatmap_DF_to_Plot["altitude (km)"].nunique() == 1):
        print()
        print("\033[1mWARNING: multiple altitudes were supplied in the input dataframe, therefore only the map for the maximum altitude will be plotted!\033[0m")
        print()
        heatmap_DF_to_Plot = heatmap_DF_to_Plot.query(f"`altitude (km)` == {heatmap_DF_to_Plot['altitude (km)'].max()}")

    if hue_range is None:
        hue_range = (0,heatmap_DF_to_Plot["edose"].max())

    ############################ creating background world map and dose image
    currentFigure = plt.gcf()

    currentFigure.set_figheight(10)
    currentFigure.set_figwidth(10)

    #heatmap_DF_to_Plot = pd.read_csv(file_path_to_read, delimiter=',')
    heatmap_DF_to_Plot["SEU (Upsets/hr/Gb)"] = heatmap_DF_to_Plot["SEU"] * (60.0 * 60.0) * 1e9
    heatmap_DF_to_Plot["SEL (Latch-ups/hr/device)"] = heatmap_DF_to_Plot["SEL"] * (60.0 * 60.0)
    if plot_longitude_east is False:
        heatmap_DF_to_Plot["longitudeTranslated"] = heatmap_DF_to_Plot["longitude"].apply(lambda x:x-360.0 if x > 180.0 else x)
    else:
        heatmap_DF_to_Plot["longitudeTranslated"] = heatmap_DF_to_Plot["longitude"]

    scatterPlotAxis = sns.scatterplot(data=heatmap_DF_to_Plot,x="longitudeTranslated",y="latitude",
                    hue=dose_type, hue_norm=hue_range, palette=palette,
                    zorder=10,
                    marker="s",s=heatmap_s,edgecolor=edgecolor,
                    legend=False,
                    )#ax=axToPlotOn)

    if plot_colorbar is True:
        norm = plt.Normalize(hue_range[0], hue_range[1])
        sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
        sm.set_array([])

        # Remove the legend and add a colorbar
        #scatterPlotAxis.get_legend().remove()
        #colorbar = scatterPlotAxis.figure.colorbar(sm,label=legend_label,shrink=0.4)
        colorbar = scatterPlotAxis.figure.colorbar(sm,label=legend_label,orientation="horizontal")
    else:
        colorbar = None

    plt.ylim([-90,90])
    plt.xlim([-175,180])
    plt.grid(True)
    plt.xlabel("Longitude (degrees)")
    plt.ylabel("Latitude (degrees)")

    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world.plot(color="None",edgecolor="black",lw=0.35,ax=scatterPlotAxis,zorder=20)
    if plot_longitude_east is True:
        world['geometry'] = world['geometry'].translate(xoff=360)
        world.plot(color="None",edgecolor="black",lw=0.35,ax=scatterPlotAxis,zorder=20)
        plt.xlim([0,355])

    ####################################################################

    #plt.legend(title=dose_type,loc="center left",bbox_to_anchor=(1.1,0.5))

    return scatterPlotAxis, colorbar

def plot_dose_map(map_to_plot,
                  plot_title=None,
                  plot_contours=True,
                  levels=3,
                    **kwargs):

    #altitude_to_plot_in_km = altitude_to_plot_in_kft * 0.3048

    axis_no_colorbar, colorbar = create_single_dose_map_plot_plt(map_to_plot,
                                                     **kwargs)

    plt.title(plot_title)

    if plot_contours is True:
        plot_dose_map_contours(map_to_plot,levels=levels,**kwargs)

    return axis_no_colorbar, colorbar

def add_colorbar_to_plot(hue_range, palette, legend_label, scatterPlotAxis=None):
    norm = plt.Normalize(hue_range[0], hue_range[1])
    sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    #scatterPlotAxis.get_legend().remove()
    if scatterPlotAxis is None:
        colorbar = plt.colorbar(sm,label=legend_label,shrink=0.4)
    else:
        colorbar = scatterPlotAxis.figure.colorbar(sm,label=legend_label,shrink=0.4)
    return colorbar

def create_single_dose_map_plotly(DF_to_use,
                              selected_altitude_in_km,
                              **kwargs):

    if selected_altitude_in_km is not None:
        DF_to_use = DF_to_use[round(DF_to_use["altitude (km)"],4) == selected_altitude_in_km]

    if len(DF_to_use) == 0:
        raise Exception("Error: specified altitude in kilometers did not match any of the altitudes in kilometers in the inputted DataFrame.")

    doseRateMap = px.scatter(DF_to_use, x="longitude",y="latitude",color="adose",
                            symbol_sequence=["square"],
                            range_y=[-90,90],
                            range_x=[0,360],
                            **kwargs)

    doseRateMap.update_traces(marker={'size': 10})
    doseRateMap.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1, range=[-90,90]))
    doseRateMap.update_xaxes(range=[0,360])
    doseRateMap.update_yaxes(range=[-90,90])
    doseRateMap.update_layout(xaxis_scaleanchor="y")

    doseRateMap.update_layout(autosize=False,
                            width=800,
                            height=600)

    doseRateMap.show()

    return doseRateMap