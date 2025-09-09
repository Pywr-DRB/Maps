"""
Contains all plotting functions used for Pywr-DRB model assessments, including:


"""
import geopandas as gpd
from shapely import ops
from shapely.geometry import Point, LineString, MultiLineString
from shapely.geometry import box

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import contextily as cx

mcm_to_mg = 264.17
mg_to_mcm = 1 / mcm_to_mg

dpi=450


spatial_data_dir = "./data/DRB_spatial/"
fig_dir = "./figures/Lin_etal_preprint/"



mcm_to_mg = 264.17
mg_to_mcm = 1 / mcm_to_mg

### Map of major DRB nodes
def make_DRB_map(fig_dir=fig_dir, 
                 scale_reservoirs_by_capacity=False,
                 use_basemap = True,
                 plot_temp_lstm_inputs=True,
                 plot_salinity_lstm_inputs=True,
                 plot_salinity_normal_range=True,
                 plot_salinity_target=True,
                 plot_tributaries = True,
                 plot_flow_requirements = True,
                 plot_lordville = True,
                 annotate_state_boundaries=True,
                 annotate_nyc_reservoirs=True,
                 annotage_drbc_lb_reservoirs=True,
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



    ### create map figures
    fig, ax = plt.subplots(1, 1, figsize=(6, 10))
    

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
    list_drbc_lb_reservoirs = [
        'reservoir_beltzvilleCombined', 
        'reservoir_blueMarsh', 
        'reservoir_fewalter'
        ]
    
    
    
    for r in reservoirs['name']:
        
        # for nyc reservoirs, plot in red
        if r in list_nyc_reservoirs:
            color = 'firebrick'
        # for drbc lb reservoirs, plot in purple
        elif r in list_drbc_lb_reservoirs:
            color = 'darkorchid'
        
        # ignore all others
        else:
            continue
                
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

        flow_reqs.plot(ax=ax, color='mediumseagreen', edgecolor='k', markersize=250, zorder=2.1, marker='*')
    
    # Lordville gauge location
    if plot_lordville:
        lordville = major_nodes.loc[major_nodes['name'] == 'link_delLordville']
        lordville = gpd.GeoDataFrame(lordville,
                                    geometry=gpd.points_from_xy(lordville.long, lordville.lat,
                                                                crs=crs_nodedata)).to_crs(crs)
        lordville.plot(ax=ax, color='gold', 
                       edgecolor='k', markersize=250, zorder=2.1, marker='*')


    # Temporary LSTM input locations
    if plot_temp_lstm_inputs:
        # two input locations are:
        # 1. cannonsville release gauge
        # 2. east branch delaware outlet
        Qin_a = major_nodes.loc[major_nodes['name'] == 'link_01425000']
        Qin_a = gpd.GeoDataFrame(Qin_a,
                                geometry=gpd.points_from_xy(Qin_a.long, Qin_a.lat,
                                                            crs=crs_nodedata)).to_crs(crs)

        # second point has coordinates:
        coords = (-75.10, 42.0)
        
        Qin_b = gpd.GeoDataFrame({'name': ['Qin_b'], 'long': [coords[0]], 'lat': [coords[1]]},
                                  geometry=gpd.points_from_xy([coords[0]], [coords[1]],
                                                              crs=crs_nodedata)).to_crs(crs)
        Qin_a.plot(ax=ax, color='blue', edgecolor='k', markersize=100, zorder=2.1, marker='D')
        Qin_b.plot(ax=ax, color='brown', edgecolor='k', markersize=100, zorder=2.1, marker='D')


    if plot_salinity_lstm_inputs:
        # two input locations are:
        # 1. trenton gage
        # 2. schuylkill outlet
        Qin_c = major_nodes.loc[major_nodes['name'] == 'link_delTrenton']
        
        Qin_c = gpd.GeoDataFrame(Qin_c,
                                geometry=gpd.points_from_xy(Qin_c.long, Qin_c.lat,
                                                            crs=crs_nodedata)).to_crs(crs)

        Qin_d = major_nodes.loc[major_nodes['name'] == 'link_outletSchuylkill']
        
        # for the schuylkill outlet, move the coordinate to the west for better visibility
        Qin_d['lat'] = float(Qin_d['lat'])
        Qin_d['long'] = float(Qin_d['long']) - 0.1
        
        Qin_d = gpd.GeoDataFrame(Qin_d,
                                geometry=gpd.points_from_xy(Qin_d.long, Qin_d.lat,
                                                            crs=crs_nodedata)).to_crs(crs)
        Qin_c.plot(ax=ax, color='mediumseagreen', edgecolor='k', markersize=100, zorder=2.1, marker='D')
        Qin_d.plot(ax=ax, color='darkgreen', edgecolor='k', markersize=100, zorder=2.1, marker='D')

    if plot_salinity_normal_range:
        # northern and southern bounds of normal salinity range
        range_lat_bounds = [39.628079, 39.878107]
        
        # convert these latitudes to the map crs
        point1 = gpd.GeoDataFrame({'name': ['point1'], 'long': [-75.477484], 'lat': [range_lat_bounds[0]]},
                                  geometry=gpd.points_from_xy([-75.477484], [range_lat_bounds[0]],
                                                              crs=crs_nodedata)).to_crs(crs)
        point2 = gpd.GeoDataFrame({'name': ['point2'], 'long': [-75.477484], 'lat': [range_lat_bounds[1]]},
                                  geometry=gpd.points_from_xy([-75.477484], [range_lat_bounds[1]],
                                                              crs=crs_nodedata)).to_crs(crs)
        
        
        # get the drb mainstem line geometry which is within these lat bounds
        # must trim the line to be within the lat bounds
        bounds = mainstem.total_bounds  # [minx, miny, maxx, maxy]
        clip_box = box(bounds[0], point1.geometry.y.iloc[0], bounds[2], point2.geometry.y.iloc[0])
        
        # Clip the geometry to the latitude range
        mainstem_range_segment = gpd.clip(mainstem, clip_box)
    
        mainstem_range_segment.plot(ax=ax, color='cyan', lw=6, alpha=0.5, zorder=1.2) 
        
    if plot_salinity_target:
        # make a single point for the target location 
        coords = (-75.477484, 39.752665)
        target_point = gpd.GeoDataFrame({'name': ['target'], 'long': [coords[0]], 'lat': [coords[1]]},
                                  geometry=gpd.points_from_xy([coords[0]], [coords[1]],
                                                              crs=crs_nodedata)).to_crs(crs)
        target_point.plot(ax=ax, color='gold', edgecolor='k', 
                          markersize=150, zorder=2.3, marker='*')
        
    ### add state boundaries
    if use_basemap:
        states.plot(ax=ax, color='none', edgecolor='0.5', lw=0.7, zorder=0)
    else:
        states.plot(ax=ax, color='0.95', edgecolor='0.5', lw=0.7, zorder=0)


    ### map limits
    ax.set_xlim([-8.517e6, -8.197e6])  # -8.215e6])
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
        
    if plot_lordville:
        plt.annotate('Lordville', xy=(-8.410e6, 5.143e6), 
                     ha='center', va='center', 
                     fontsize=fontsize, color='gold', fontweight='bold')

    if plot_flow_requirements:
        fontcolor = 'mediumseagreen'
        plt.annotate('Montague', xy=(-8.284e6, 5.055e6), 
                     ha='center', va='center', fontsize=fontsize, color=fontcolor,
                    fontweight='bold')
        plt.annotate('Trenton', xy=(-8.353e6, 4.894e6), 
                     ha='center', va='center', fontsize=fontsize, color=fontcolor,
                    fontweight='bold')





    ### legend
    axin = ax.inset_axes([0.58, 0.006, 0.41, 0.25])
    axin.set_xlim([0, 1])
    axin.set_ylim([0, 1])
    ### mainstem
    axin.plot([0.05, 0.15], [0.93, 0.93], color='navy', lw=3)
    axin.annotate('Delaware River', xy=(0.18, 0.93), ha='left', va='center', color='k', fontsize=fontsize)

    # tributaries
    if plot_tributaries:
        axin.plot([0.05, 0.15], [0.83, 0.83], color='cornflowerblue', lw=2)
        axin.annotate('Tributary', xy=(0.18, 0.83), ha='left', va='center', color='k', fontsize=fontsize)

    # DRB boundary
    axin.plot([0.05, 0.15], [0.73, 0.73], color='k', lw=1)
    axin.annotate('Basin Boundary', xy=(0.18, 0.73), ha='left', va='center', color='k', fontsize=fontsize)

    ### Minimum flow targets

    if plot_flow_requirements:
        axin.scatter([0.1], [0.53], color='mediumseagreen', edgecolor='k', s=200, marker='*')
        axin.annotate('Flow Target', xy=(0.18, 0.53), ha='left', va='center', color='mediumseagreen', fontweight='bold',
                    fontsize=fontsize)

    ### NYC Reservoirs
    axin.scatter([0.1], [0.43], color='firebrick', edgecolor='k', s=100)
    axin.annotate('NYC Reservoir', xy=(0.18, 0.43), ha='left', va='center', color='firebrick', fontweight='bold',
                  fontsize=fontsize)
    ### Non-NYC Reservoirs
    axin.scatter([0.1], [0.33], color='sandybrown', edgecolor='k', s=100)
    axin.annotate('Non-NYC Reservoir', xy=(0.18, 0.33), ha='left', va='center', color='k', fontsize=fontsize)
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

    # ### basemap - this is slow and breaks sometimes, if so just try later
    if use_basemap:
        cx.add_basemap(ax=ax, alpha=0.5, attribution_size=6, 
                       source=cx.providers.CartoDB.Positron)
        figname = f'{fig_dir}/static_map_withbasemap_positron.png'
    else:
        figname = f'{fig_dir}static_map.png'

    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig(figname, bbox_inches='tight', dpi=dpi)


if __name__ == "__main__":
    make_DRB_map(plot_flow_requirements=False, 
                 plot_salinity_lstm_inputs=True,
                 plot_temp_lstm_inputs=True,
                 scale_reservoirs_by_capacity=False,)