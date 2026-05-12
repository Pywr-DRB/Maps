"""
Contains all plotting functions used for Pywr-DRB model assessments, including:


"""
import os
import geopandas as gpd
from shapely import ops
from shapely.geometry import Point, LineString, MultiLineString
from shapely.geometry import box
from matplotlib.markers import MarkerStyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import contextily as cx
from pyproj import Transformer

mcm_to_mg = 264.17
mg_to_mcm = 1 / mcm_to_mg

dpi=450


spatial_data_dir = "./data/DRB_spatial/"
fig_dir = "./figures/Amestoy_etal_2026/"
os.makedirs(fig_dir, exist_ok=True)


### Map of major DRB nodes
def make_DRB_map(fig_dir=fig_dir,
                 scale_reservoirs_by_capacity=False,
                 use_basemap = True,
                 basemap_source = cx.providers.Esri.WorldShadedRelief,
                 plot_tributaries = True,
                 plot_flow_requirements = True,
                 annotate_state_boundaries=True,
                 annotate_nyc_reservoirs=True,
                 drb_highlight=None,
                 units='MG'):

    ### set crs consistent with contextily basemap
    crs = 'EPSG:3857'
    crs_nodedata = 4386
    crs_longlat = 'EPSG:4326'


    # # Reservoir data
    reservoir_data = pd.read_csv('./data/istarf_conus.csv', sep=',')

    ### Load general shapefiles
    drb_boundary = gpd.read_file(f'{spatial_data_dir}/DRB_shapefiles/drb_bnd_polygon.shp').to_crs(crs)
    states = gpd.read_file(f'{spatial_data_dir}/states/tl_2010_us_state10.shp').to_crs(crs)
    nhd = gpd.read_file(f'{spatial_data_dir}/NHD_0204/Shape/NHDFlowline.shp').to_crs(crs)
    drc = gpd.read_file(f'{spatial_data_dir}/NHD_NJ/DRCanal.shp').to_crs(crs)
    nyc = gpd.read_file('./data/nybb_26a/nybb.shp').to_crs(crs).dissolve()

    ### load drb node info into geodataframes
    major_nodes = gpd.read_file(f'{spatial_data_dir}/model_components/drb_model_major_nodes.csv', sep=',')
    reservoirs = major_nodes.loc[major_nodes['type'] == 'reservoir']
    reservoirs = gpd.GeoDataFrame(reservoirs,
                                  geometry=gpd.points_from_xy(reservoirs.long, reservoirs.lat,
                                                              crs=crs_nodedata)).to_crs(crs)

    ### get river network from NHD
    mainstem = nhd.loc[nhd['gnis_name'] == 'Delaware River']
    ## for river/stream objects, merge rows into single geometry to avoid white space on plot
    multi_linestring = MultiLineString([ls for ls in mainstem['geometry'].values])
    merged_linestring = ops.linemerge(multi_linestring)
    mainstem = gpd.GeoDataFrame({'geometry': [merged_linestring]})

    # smooth the line
    mainstem['geometry'] = mainstem['geometry'].apply(lambda x: x.simplify(tolerance=100))

    ### rough lines for NYC Delaware Aqueduct system
    ### cannonsville to rondout
    lines = [LineString([Point(-75.37462, 42.065872), Point(-74.4296, 41.79926)])]
    ### pepacton to rondout
    lines.append(LineString([Point(-74.965531, 42.073603), Point(-74.4296, 41.79926)]))
    ### neversink to rondout
    lines.append(LineString([Point(-74.643266, 41.821286), Point(-74.4296, 41.79926)]))
    ### rondout to west branch reservoir
    lines.append(LineString([Point(-74.4296, 41.79926), Point(-73.69541, 41.41176)]))
    ### west branch reservoir to kensico reservoir
    lines.append(LineString([Point(-73.69541, 41.41176), Point(-73.7659656, 41.0737078)]))
    ### kensico reservoir to hillside reservoir
    lines.append(LineString([Point(-73.7659656, 41.0737078), Point(-73.8693806, 40.90715556)]))

    ### convert projection
    crs_longlat = 'EPSG:4326'
    aqueducts = gpd.GeoDataFrame({'geometry': lines}, crs=crs_longlat)
    aqueducts = aqueducts.to_crs(crs)

    ### create map figures
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    

    ### highlight DRB basin (above basemap, below map features)
    if drb_highlight == 'inside':
        drb_boundary.plot(ax=ax, color='khaki', edgecolor='none',
                          alpha=0.35, zorder=0.2)
    elif drb_highlight == 'outside':
        _y_range = 5.235e6 - 4.75e6
        _x_center = (-8.517e6 + -8.197e6) / 2
        _x_shift = 0.035e6
        mask_bbox = box(_x_center - _y_range/2 + _x_shift,
                        4.75e6,
                        _x_center + _y_range/2 + _x_shift,
                        5.235e6)
        outside_geom = mask_bbox.difference(drb_boundary.unary_union)
        gpd.GeoDataFrame({'geometry': [outside_geom]}, crs=crs).plot(
            ax=ax, color='0.85', edgecolor='none', alpha=0.6, zorder=0.2)

    ### plot drb boundary
    drb_boundary.plot(ax=ax, color='none', edgecolor='k', lw=0.5, zorder=0.9)

    ### plot river network
    mainstem.plot(ax=ax, color='navy', lw=2, zorder=1.1)

    # Plot tributaries
    if plot_tributaries:
        
        # get all other tributary streams containing a Pywr-DRB reservoir, 
        # or downstream of one. Note that 2 different regions have Tulpehocken Creek - use only the correct one.
        trib_names = ['West Branch Delaware River', 'East Branch Delaware River', 'Neversink River',
                    'Mongaup River', 'Lackawaxen River', 'West Branch Lackawaxen River', 'Wallenpaupack Creek',
                    'Lehigh River', 'Shohola Creek', 'Pohopoco Creek', 'Merrill Creek', 'Musconetcong River',
                    'Pohatcong Creek', 'Tohickon Creek', 'Assunpink Creek', 'Schuylkill River', 'Maiden Creek',
                    'Tulpehocken Creek', 'Still Creek', 'Little Schuylkill River',
                    'Perkiomen Creek']

        tribs = []
        for trib_name in trib_names:
            trib = nhd.loc[[(n == trib_name) and ((n != 'Tulpehocken Creek') or ('02040203' in c)) for n, c in
                            zip(nhd['gnis_name'], nhd['reachcode'])]]
            multi_linestring = MultiLineString([ls for ls in trib['geometry'].values])
            merged_linestring = ops.linemerge(multi_linestring)
            trib = gpd.GeoDataFrame({'geometry': [merged_linestring], 'name': trib_name})

            # smooth the line
            trib['geometry'] = trib['geometry'].apply(lambda x: x.simplify(tolerance=100))

            tribs.append(trib)

        for trib in tribs:
            trib.plot(ax=ax, color='cornflowerblue', lw=1, zorder=1)

    # reservoirs
    list_nyc_reservoirs = ('reservoir_cannonsville', 'reservoir_pepacton', 'reservoir_neversink')

    for r in reservoirs['name']:
        color = 'firebrick' if r in list_nyc_reservoirs else 'sandybrown'

        if scale_reservoirs_by_capacity:
            r_abbrev = r.split('_')[1]
            if units == 'MG':
                try:
                    s = 50 + reservoir_data['Adjusted_CAP_MG'].loc[reservoir_data['reservoir'] == r_abbrev].iloc[0] / 1000 * 2
                except:
                    s = 50
            elif units == 'MCM':
                try:
                    s = 50 + reservoir_data['Adjusted_CAP_MG'].loc[reservoir_data['reservoir'] == r_abbrev].iloc[0] * mg_to_mcm / 2
                except:
                    s = 50
        else:
            s = 100

        reservoirs.loc[reservoirs['name'] == r].plot(ax=ax, color=color, edgecolor='k', markersize=s, zorder=2)

    # Montague and Trenton flow requirement locations
    if plot_flow_requirements:

        flow_reqs = major_nodes.loc[major_nodes['type'] == 'regulatory']
        flow_reqs = gpd.GeoDataFrame(flow_reqs,
                                    geometry=gpd.points_from_xy(flow_reqs.long, flow_reqs.lat,
                                                                crs=crs_nodedata)).to_crs(crs)

        montague_color = 'lightgreen'
        trenton_color = 'darkgreen'
        flow_reqs.loc[flow_reqs['name'] == 'link_delMontague'].plot(
            ax=ax, color=montague_color, edgecolor='k', markersize=250, zorder=2.1, marker='*')
        flow_reqs.loc[flow_reqs['name'] == 'link_delTrenton'].plot(
            ax=ax, color=trenton_color, edgecolor='k', markersize=250, zorder=2.1, marker='*')

    ### NYC tunnel system
    aqueducts.plot(ax=ax, color='darkmagenta', lw=2, zorder=1.2, ls=':')

    ### plot NJ diversion - Delaware & Raritan Canal
    drc.plot(ax=ax, color='darkmagenta', lw=2, zorder=1.2, ls=':')

    ### add state boundaries
    if use_basemap:
        states.plot(ax=ax, color='none', edgecolor='0.5', lw=0.7, zorder=0)
    else:
        states.plot(ax=ax, color='0.95', edgecolor='0.5', lw=0.7, zorder=0)

    ### highlight NYC
    nyc.plot(ax=ax, color='goldenrod', edgecolor='goldenrod', lw=0.6, alpha=0.75, zorder=1.9)
    nyc_centroid = nyc.geometry.iloc[0].centroid
    plt.annotate('NYC', xy=(nyc_centroid.x + 0.008e6, nyc_centroid.y + 0.004e6),
                 ha='center', va='center', fontsize=10,
                 color='k', fontweight='bold', zorder=2.2)


    ### map limits
    # Make square by adjusting x limits to match y range
    y_range = 5.235e6 - 4.75e6  # 0.485e6
    x_center = (-8.517e6 + -8.197e6) / 2  # -8.357e6
    # Shift view east (increase x limits to show more eastern content)
    x_shift = 0.035e6
    ax.set_xlim([x_center - y_range/2 + x_shift, x_center + y_range/2 + x_shift])
    ax.set_ylim([4.75e6, 5.235e6])

    ### annotations
    if annotate_state_boundaries:
        fontsize = 10
        fontcolor = '0.5'
        plt.annotate('New York', xy=(-8.46e6, 5.168e6), ha='center', va='center', fontsize=fontsize, color=fontcolor)
        plt.annotate('Pennsylvania', xy=(-8.46e6, 5.152e6), ha='center', va='center', fontsize=fontsize, color=fontcolor)
        plt.annotate('New York', xy=(-8.245e6, 5.032e6), rotation=-31, ha='center', va='center', fontsize=fontsize,
                    color=fontcolor)
        plt.annotate('New Jersey', xy=(-8.257e6, 5.019e6), rotation=-31, ha='center', va='center', fontsize=fontsize,
                    color=fontcolor)
        plt.annotate('Pennsylvania', xy=(-8.48e6, 4.833e6), ha='center', va='center', fontsize=fontsize, color=fontcolor)
        plt.annotate('Maryland', xy=(-8.48e6, 4.817e6), ha='center', va='center', fontsize=fontsize, color=fontcolor)
        plt.annotate('Delaware', xy=(-8.43e6, 4.8e6), rotation=-85, ha='center', va='center', fontsize=fontsize,
                    color=fontcolor)
        plt.annotate('Pennsylvania', xy=(-8.359e6, 5.03e6), rotation=50, ha='center', va='center', fontsize=fontsize,
                    color=fontcolor)
        plt.annotate('New Jersey', xy=(-8.337e6, 5.024e6), rotation=50, ha='center', va='center', fontsize=fontsize,
                    color=fontcolor)

    if annotate_nyc_reservoirs:
        fontcolor = 'firebrick'
        plt.annotate('Cannonsville', xy=(-8.404e6, 5.1835e6), ha='center', va='center', fontsize=fontsize, color=fontcolor,
                    fontweight='bold')
        plt.annotate('Pepacton', xy=(-8.306e6, 5.168e6), ha='center', va='center', fontsize=fontsize, color=fontcolor,
                    fontweight='bold')
        plt.annotate('Neversink', xy=(-8.268e6, 5.143e6), 
                     ha='center', va='center', 
                     fontsize=fontsize, color=fontcolor, fontweight='bold')

    if plot_flow_requirements:
        montague_label = plt.annotate('Montague', xy=(-8.284e6, 5.055e6),
                     ha='center', va='center', fontsize=fontsize, color='lightgreen',
                    fontweight='bold')
        montague_label.set_path_effects([
            path_effects.Stroke(linewidth=1.2, foreground='black'),
            path_effects.Normal()])
        plt.annotate('Trenton', xy=(-8.353e6, 4.894e6),
                     ha='center', va='center', fontsize=fontsize, color='darkgreen',
                    fontweight='bold')

    # Diversion labels
    fontcolor = 'darkmagenta'
    plt.annotate('NYC\nDiversion', xy=(-8.25e6, 5.085e6), ha='center', va='center',
                 fontsize=fontsize, color=fontcolor, fontweight='bold')
    plt.annotate('NJ\nDiversion', xy=(-8.315e6, 4.959e6), ha='center', va='center',
                 fontsize=fontsize, color=fontcolor, fontweight='bold')





    ### Custom legend
    # Determine how many items will be included
    n_legend_items = 7 # mainstem, boundary, nyc res, non-nyc res, diversions, reservoir cap scale, NYC

    n_legend_items += int(plot_tributaries)
    n_legend_items += 2 * int(plot_flow_requirements)
    
    # current y position for legend items
    yi = 0.95
    item_spacing = 1/n_legend_items
    fontsize = 9

    # make inset axes for legend
    axin_height = 0.03 * n_legend_items + 0.06
    axin_width = 0.35
    axin_x0 = 0.65
    axin_y0 = 0.0
    axin = ax.inset_axes([axin_x0, axin_y0, 
                          axin_width, axin_height])  # [x0, y0, width, height]
    axin.set_xlim([0, 1])
    axin.set_ylim([0, 1])
    
    
    ### mainstem
    axin.plot([0.05, 0.15], [yi, yi], color='navy', lw=3)
    axin.annotate('Delaware River', xy=(0.18, yi), ha='left', va='center', color='k', fontsize=fontsize)

    # tributaries
    if plot_tributaries:
        yi -= item_spacing
        axin.plot([0.05, 0.15], [yi, yi], color='cornflowerblue', lw=2)
        axin.annotate('Tributary', xy=(0.18, yi), ha='left', va='center', color='k', fontsize=fontsize)

    # DRB boundary
    yi -= item_spacing
    axin.plot([0.05, 0.15], [yi, yi], color='k', lw=1)
    axin.annotate('Basin Boundary', xy=(0.18, yi), ha='left', va='center', color='k', fontsize=fontsize)

    ### Diversions
    yi -= item_spacing
    axin.plot([0.05, 0.15], [yi, yi], color='darkmagenta', lw=2, ls=':')
    axin.annotate('Interbasin Transfer', xy=(0.18, yi), ha='left', va='center', color='k', fontsize=fontsize)

    ### NYC Reservoirs
    yi -= item_spacing
    axin.scatter([0.1], [yi], color='firebrick', edgecolor='k', s=100)
    axin.annotate('NYC Reservoir', xy=(0.18, yi), 
                  ha='left', va='center', color='k', 
                  fontsize=fontsize)

    # Non-NYC Reservoirs
    yi -= item_spacing
    axin.scatter([0.1], [yi], color='sandybrown', edgecolor='k', s=100)
    axin.annotate('Non-NYC Reservoir', xy=(0.18, yi),
                  ha='left', va='center', color='k',
                  fontsize=fontsize)

    ### NYC boundary
    yi -= item_spacing
    axin.add_patch(mpatches.Rectangle((0.05, yi - 0.02), 0.10, 0.04,
                                      facecolor='goldenrod', edgecolor='goldenrod',
                                      lw=0.6, alpha=0.75))
    axin.annotate('New York City', xy=(0.18, yi),
                  ha='left', va='center', color='k', fontsize=fontsize)

    ### Minimum flow targets
    if plot_flow_requirements:
        yi -= item_spacing
        axin.scatter([0.1], [yi], color='lightgreen', edgecolor='k', s=200, marker='*')
        axin.annotate('Montague Flow Target', xy=(0.18, yi),
                      ha='left', va='center', color='k',
                      fontsize=fontsize)
        yi -= item_spacing
        axin.scatter([0.1], [yi], color='darkgreen', edgecolor='k', s=200, marker='*')
        axin.annotate('Trenton Flow Target', xy=(0.18, yi),
                      ha='left', va='center', color='k',
                      fontsize=fontsize)

    ### Non-NYC Reservoirs
    # axin.scatter([0.1], [0.33], color='sandybrown', edgecolor='k', s=100)
    # axin.annotate('Non-NYC Reservoir', xy=(0.18, 0.33), ha='left', va='center', color='k', fontsize=fontsize)


    ### marker size for reservoirs
    # axin.annotate('Reservoir Capacity', xy=(0.05, 0.3), ha='left', va='center', color='k', fontsize=fontsize)    
    
    if scale_reservoirs_by_capacity:
        axin.scatter([0.15], [0.18], color='0.5', edgecolor='k', s=50 + 10  / 2)
        axin.scatter([0.45], [0.18], color='0.5', edgecolor='k', s=50 + 200  / 2)
        axin.scatter([0.8], [0.18], color='0.5', edgecolor='k', s=50 + 530  / 2)
        axin.annotate('10', xy=(0.15, 0.05), ha='center', va='center', color='k', fontsize=fontsize)
        axin.annotate('200', xy=(0.45, 0.05), ha='center', va='center', color='k', fontsize=fontsize)
        axin.annotate('530 MCM', xy=(0.8, 0.05), ha='center', va='center', color='k', fontsize=fontsize)
    ### clean up
    axin.set_xticks([])
    axin.set_yticks([])
    axin.patch.set_alpha(0.9)

    ### Long/lat tick labels.
    # Map axes are in EPSG:3857 (Web Mercator). Pick lat/long values in the
    # visible extent, project them to Web Mercator with pyproj for accurate
    # placement, then label the ticks in degrees.
    to_webmercator = Transformer.from_crs(crs_longlat, crs, always_xy=True)
    lon_ticks_deg = [-76, -75, -74, -73]
    lat_ticks_deg = [40, 41, 42]
    x_ticks = [to_webmercator.transform(lon, 40.0)[0] for lon in lon_ticks_deg]
    y_ticks = [to_webmercator.transform(-75.0, lat)[1] for lat in lat_ticks_deg]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{abs(lon)}°W' for lon in lon_ticks_deg], fontsize=9)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{lat}°N' for lat in lat_ticks_deg], fontsize=9)
    ax.tick_params(axis='both', which='major', direction='out', length=4, pad=2)
    
    # ### basemap - this is slow and breaks sometimes, if so just try later
    if use_basemap:
        cx.add_basemap(ax=ax, alpha=0.75, attribution_size=6,
                       source=basemap_source)

        # Add north arrow in upper left corner
        arrow_x = -8.55e6  # Shifted further left
        arrow_y = 5.17e6   # Shifted down from 5.20e6
        arrow_length = 0.025e6  # Reduced from 0.04e6 to ensure top is visible

        arrow = mpatches.FancyArrow(arrow_x, arrow_y, 0, arrow_length,
                                    width=arrow_length*0.15,
                                    head_width=arrow_length*0.4,
                                    head_length=arrow_length*0.25,
                                    fc='black', ec='black',
                                    linewidth=1, zorder=10)
        ax.add_patch(arrow)

        ax.text(arrow_x, arrow_y + arrow_length + 0.01e6, 'N',
                fontsize=12, fontweight='bold', ha='center', va='bottom',
                zorder=10)

        # Add scale bar in bottom left corner
        scale_x = -8.54e6  # Shifted further left
        scale_y = 4.78e6   # Near bottom edge

        # Calculate scale bar length accounting for Web Mercator distortion
        lat_center = 40.5  # degrees (map center latitude)
        scale_factor = 1 / np.cos(np.radians(lat_center))
        scale_length_km = 50  # kilometers
        scale_length_map = scale_length_km * 1000 * scale_factor

        # Draw horizontal baseline with tick marks at 0, 25, 50 km
        n_ticks = 2
        tick_height = 0.006e6
        ax.plot([scale_x, scale_x + scale_length_map], [scale_y, scale_y],
                color='black', lw=1.5, solid_capstyle='butt', zorder=10)
        for i in range(n_ticks + 1):
            tx = scale_x + i * (scale_length_map / n_ticks)
            ax.plot([tx, tx], [scale_y, scale_y + tick_height],
                    color='black', lw=1.5, zorder=10)
            tick_km = int(round(scale_length_km * i / n_ticks))
            label = f'{tick_km} km' if i == n_ticks else f'{tick_km}'
            ax.text(tx, scale_y - 0.004e6, label,
                    fontsize=8, ha='center', va='top', zorder=10)

        # Add inlay map showing regional context in upper right
        axin_inlay = ax.inset_axes([0.70, 0.75, 0.28, 0.23])  # [x0, y0, width, height]

        # Filter states to NY, PA, NJ, DE, MD
        state_abbrevs = ['NY', 'PA', 'NJ', 'DE', 'MD']
        try:
            regional_states = states[states['STUSPS10'].isin(state_abbrevs)]
        except KeyError:
            # Fallback if column name is different
            print("Warning: STUSPS10 column not found, trying alternatives")
            for col in ['STUSPS', 'STATE_ABBR', 'USPS']:
                if col in states.columns:
                    regional_states = states[states[col].isin(state_abbrevs)]
                    break

        # Plot regional states
        regional_states.plot(ax=axin_inlay, color='0.9', edgecolor='black',
                             linewidth=0.8, zorder=1)

        # Highlight DRB basin
        drb_boundary.plot(ax=axin_inlay, color='steelblue', alpha=0.6,
                          edgecolor='navy', linewidth=1.5, zorder=2)

        # Set extent with padding
        bounds = regional_states.total_bounds  # [minx, miny, maxx, maxy]
        padding = 0.05 * (bounds[2] - bounds[0])
        axin_inlay.set_xlim([bounds[0] - padding, bounds[2] + padding])
        axin_inlay.set_ylim([bounds[1] - padding, bounds[3] + padding])

        # Clean up axes
        axin_inlay.set_xticks([])
        axin_inlay.set_yticks([])
        axin_inlay.set_aspect('equal')
        for spine in axin_inlay.spines.values():
            spine.set_linewidth(1.5)

        # convert basemap soure to str
        basemap_source_str = basemap_source.name
        highlight_suffix = f'_drb{drb_highlight}' if drb_highlight else ''
        figname = f'{fig_dir}/static_map_withbasemap_{basemap_source_str}{highlight_suffix}.png'
        plt.savefig(figname, bbox_inches='tight', dpi=dpi)

    else:
        figname = f'{fig_dir}static_map.png'


if __name__ == "__main__":
    # Generate map with default Esri.WorldShadedRelief basemap
    make_DRB_map(plot_flow_requirements=True,
                 scale_reservoirs_by_capacity=False,
                 drb_highlight='outside')