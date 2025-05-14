import ee
import geemap
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import os
import streamlit as st
from datetime import datetime
import tempfile

# Set page config
st.set_page_config(layout="wide", page_title="Farm Geospatial Analysis")

# App title and description
st.title("Farm Geospatial Analysis Dashboard")
st.markdown("""
This app provides geospatial analysis for farm properties using Google Earth Engine data.
Upload your farm shapefile (zipped) or use coordinates to analyze various environmental indicators.
""")

# Initialize GEE
try:
    ee.Initialize()
except Exception as e:
    st.error(f"Error initializing Google Earth Engine: {e}")
    st.info("If you're running this locally, make sure you've authenticated with GEE using 'earthengine authenticate'.")
    st.stop()

# Sidebar for inputs
st.sidebar.header("Input Options")
input_method = st.sidebar.radio("Select input method:", ["Upload Shapefile", "Enter Coordinates"])

farm_ee = None
farm_gdf = None

if input_method == "Upload Shapefile":
    uploaded_file = st.sidebar.file_uploader("Upload zipped shapefile", type=['zip'])
    
    if uploaded_file is not None:
        # Create a temporary directory to extract the zip files
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save the uploaded zip file
            zip_path = os.path.join(temp_dir, "farm_shapefile.zip")
            with open(zip_path, "wb") as f:
                f.write(uploaded_file.getvalue())
            
            # Use GeoPandas to read the shapefile from the zip
            try:
                farm_gdf = gpd.read_file(f"zip://{zip_path}")
                farm_gdf = farm_gdf.to_crs(epsg=4326)
                farm_ee = geemap.gdf_to_ee(farm_gdf)
                
                st.sidebar.success("Shapefile loaded successfully!")
                st.sidebar.map(farm_gdf)
            except Exception as e:
                st.sidebar.error(f"Error reading shapefile: {e}")
                st.sidebar.info("Make sure your zip file contains all necessary shapefile components (.shp, .shx, .dbf, etc.)")

else:  # Coordinate input
    st.sidebar.subheader("Enter Farm Polygon Coordinates")
    st.sidebar.markdown("Enter coordinates in decimal degrees (longitude, latitude)")
    
    # Create coordinate input fields
    coords = []
    col1, col2 = st.sidebar.columns(2)
    
    with col1:
        st.write("Longitude")
        lon1 = st.number_input("Point 1 Lon", value=-2.5, format="%.5f")
        lon2 = st.number_input("Point 2 Lon", value=-2.5, format="%.5f")
        lon3 = st.number_input("Point 3 Lon", value=-2.49, format="%.5f")
        lon4 = st.number_input("Point 4 Lon", value=-2.49, format="%.5f")
    
    with col2:
        st.write("Latitude")
        lat1 = st.number_input("Point 1 Lat", value=36.85, format="%.5f")
        lat2 = st.number_input("Point 2 Lat", value=36.86, format="%.5f")
        lat3 = st.number_input("Point 3 Lat", value=36.86, format="%.5f")
        lat4 = st.number_input("Point 4 Lat", value=36.85, format="%.5f")
    
    farm_coords = [
        [lon1, lat1], [lon2, lat2], [lon3, lat3], [lon4, lat4], [lon1, lat1]
    ]
    
    farm_ee = ee.Geometry.Polygon(farm_coords)
    
    # Create a GeoDataFrame for mapping
    polygon = gpd.GeoSeries([gpd.Polygon(farm_coords)], crs="EPSG:4326")
    farm_gdf = gpd.GeoDataFrame(geometry=polygon)
    
    # Show coordinates on mini-map
    st.sidebar.map(farm_gdf)

# Analysis parameters
st.sidebar.header("Analysis Parameters")
start_date = st.sidebar.date_input("Start Date", datetime(2024, 1, 1))
end_date = st.sidebar.date_input("End Date", datetime(2024, 12, 31))
cloud_threshold = st.sidebar.slider("Cloud Coverage Threshold (%)", 0, 100, 10)

# Proceed only if we have a farm boundary
if farm_ee:
    # Create tabs for different analyses
    tab1, tab2, tab3, tab4 = st.tabs(["Vegetation & Water", "Terrain & Soil", "Climate", "Conservation"])
    
    # Use a spinner to show processing
    with st.spinner("Processing data in Google Earth Engine..."):
        # Sentinel-2 data
        s2 = ee.ImageCollection('COPERNICUS/S2_SR') \
            .filterBounds(farm_ee) \
            .filterDate(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_threshold)) \
            .median()
        
        # Calculate indices
        ndvi = s2.normalizedDifference(['B8', 'B4']).rename('NDVI')
        ndvi_vis = {'min': 0.1, 'max': 0.5, 'palette': ['white', 'green']}
        
        ndwi = s2.normalizedDifference(['B3', 'B8']).rename('NDWI')
        ndwi_vis = {'min': -1, 'max': 1, 'palette': ['white', 'blue']}
        
        # Land Cover data
        landcover = ee.Image('ESA/WorldCover/v100/2020').select('Map')
        landcover_vis = {
            'min': 10, 'max': 100,
            'palette': [
                '006400',  # Tree cover
                'ffbb22',  # Shrubland
                'ffff4c',  # Grassland
                'f096ff',  # Cropland
                'fa0000',  # Built-up
                'b4b4b4',  # Bare/sparse
                'f0f0f0',  # Snow/ice
                '0064c8',  # Water bodies
                '0096a0',  # Herbaceous wetland
                '00cf75',  # Mangroves
                'fae6a0'   # Moss/lichen
            ]
        }
        
        # Terrain data
        elevation = ee.Image('USGS/SRTMGL1_003')
        slope = ee.Terrain.slope(elevation)
        slope_vis = {'min': 0, 'max': 30, 'palette': ['white', 'black']}
        
        # Soil data
        soil = ee.Image('OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02').select('b0').clip(farm_ee)
        soil_vis = {'min': 2, 'max': 5, 'palette': ['ffffa0', 'f7fcb9', 'd9f0a3', 'addd8e', '78c679', '41ab5d', '238443', '005b29', '004b29', '012b13', '00120b']}
        
        # Climate data
        precip = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
            .filterBounds(farm_ee) \
            .filterDate(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')) \
            .sum()
        precip_vis = {'min': 0, 'max': 1000, 'palette': ['#FFFF00', 'blue']}
        
        drought = ee.ImageCollection('CSIC/SPEI/2_8') \
            .filterBounds(farm_ee) \
            .filterDate(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')) \
            .mean()
        drought_vis = {'min': -3, 'max': 2, 'palette': ['ff4819', '0095ff']}
        
        # Conservation data
        wdpa = ee.FeatureCollection('WCMC/WDPA/current/polygons')
        wdpa_vis = {'color': 'purple', 'fillColor': 'purple', 'fillOpacity': 0.3}
        
        hotspots = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017') \
            .filter(ee.Filter.eq('BIOME_NAME', 'Mediterranean Forests, Woodlands & Scrub')) \
            .filterBounds(farm_ee)
        hotspots_vis = {'color': 'red', 'fillColor': 'red', 'fillOpacity': 0.2}
        
        # NDVI time series
        ndvi_time_series = ee.ImageCollection('COPERNICUS/S2_SR') \
            .filterBounds(farm_ee) \
            .filterDate(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d')) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_threshold)) \
            .map(lambda image: image.normalizedDifference(['B8', 'B4']).rename('NDVI') \
                 .set('system:time_start', image.get('system:time_start')))
        
        def get_ndvi_stats(image):
            stats = image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=farm_ee,
                scale=10
            )
            return ee.Feature(None, {
                'date': ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
                'ndvi': stats.get('NDVI')
            })
        
        ndvi_stats = ndvi_time_series.map(get_ndvi_stats).getInfo()
        
        # Conservation overlap checks
        hotspot_overlap = hotspots.filterBounds(farm_ee).size().getInfo()
        wdpa_overlap = wdpa.filterBounds(farm_ee).size().getInfo()

    # Tab 1: Vegetation & Water
    with tab1:
        st.header("Vegetation & Water Analysis")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # NDVI Map
            m_ndvi = geemap.Map()
            m_ndvi.centerObject(farm_ee, 14)
            m_ndvi.addLayer(ndvi.clip(farm_ee), ndvi_vis, 'NDVI')
            m_ndvi.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            ndvi_mapbox = m_ndvi.to_streamlit(height=400)
            st.write("NDVI (Normalized Difference Vegetation Index)")
        
        with col2:
            # NDWI Map
            m_ndwi = geemap.Map()
            m_ndwi.centerObject(farm_ee, 14)
            m_ndwi.addLayer(ndwi.clip(farm_ee), ndwi_vis, 'NDWI')
            m_ndwi.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            ndwi_mapbox = m_ndwi.to_streamlit(height=400)
            st.write("NDWI (Normalized Difference Water Index)")
        
        # Land Cover
        st.subheader("Land Cover")
        m_lc = geemap.Map()
        m_lc.centerObject(farm_ee, 14)
        m_lc.addLayer(landcover.clip(farm_ee), landcover_vis, 'Land Cover')
        m_lc.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
        lc_mapbox = m_lc.to_streamlit(height=400)
        
        # NDVI Time Series
        st.subheader("NDVI Time Series")
        if ndvi_stats['features']:
            dates = [item['properties']['date'] for item in ndvi_stats['features'] if 'date' in item['properties'] and 'ndvi' in item['properties'] and item['properties']['ndvi'] is not None]
            ndvi_values = [item['properties']['ndvi'] for item in ndvi_stats['features'] if 'date' in item['properties'] and 'ndvi' in item['properties'] and item['properties']['ndvi'] is not None]
            
            if dates and ndvi_values:
                df = pd.DataFrame({'Date': dates, 'NDVI': ndvi_values})
                st.line_chart(df.set_index('Date'))
            else:
                st.warning("No valid NDVI data available for the selected time period.")
        else:
            st.warning("No NDVI time series data available for the selected time period and cloud coverage threshold.")
    
    # Tab 2: Terrain & Soil
    with tab2:
        st.header("Terrain & Soil Analysis")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Slope Map
            m_slope = geemap.Map()
            m_slope.centerObject(farm_ee, 14)
            m_slope.addLayer(slope.clip(farm_ee), slope_vis, 'Slope')
            m_slope.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            slope_mapbox = m_slope.to_streamlit(height=400)
            st.write("Slope (degrees)")
        
        with col2:
            # Soil Map
            m_soil = geemap.Map()
            m_soil.centerObject(farm_ee, 14)
            m_soil.addLayer(soil, soil_vis, 'Soil Organic Carbon')
            m_soil.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            soil_mapbox = m_soil.to_streamlit(height=400)
            st.write("Soil Organic Carbon")
    
    # Tab 3: Climate
    with tab3:
        st.header("Climate Analysis")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Precipitation Map
            m_precip = geemap.Map()
            m_precip.centerObject(farm_ee, 14)
            m_precip.addLayer(precip.clip(farm_ee), precip_vis, 'Precipitation')
            m_precip.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            precip_mapbox = m_precip.to_streamlit(height=400)
            st.write("Precipitation (mm)")
        
        with col2:
            # Drought Map
            m_drought = geemap.Map()
            m_drought.centerObject(farm_ee, 14)
            m_drought.addLayer(drought.clip(farm_ee), drought_vis, 'Drought (SPEI)')
            m_drought.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
            drought_mapbox = m_drought.to_streamlit(height=400)
            st.write("Drought Index (SPEI)")
    
    # Tab 4: Conservation
    with tab4:
        st.header("Conservation Analysis")
        
        # Protected Areas and Ecoregions Map
        m_conserv = geemap.Map()
        m_conserv.centerObject(farm_ee, 10)
        m_conserv.addLayer(hotspots, hotspots_vis, 'Mediterranean Ecoregions')
        m_conserv.addLayer(wdpa.filterBounds(farm_ee.buffer(10000)), wdpa_vis, 'Protected Areas')
        m_conserv.addLayer(farm_ee, {'color': 'black'}, 'Farm Boundary')
        conserv_mapbox = m_conserv.to_streamlit(height=500)
        
        # Conservation overlap info
        st.subheader("Conservation Status")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Mediterranean Ecoregion Overlap", 
                     "Within Ecoregion" if hotspot_overlap > 0 else "Not in Ecoregion")
        
        with col2:
            st.metric("Protected Area Overlap", 
                     "Within Protected Area" if wdpa_overlap > 0 else "Not in Protected Area")
        
        # Conservation recommendations based on overlap
        st.subheader("Conservation Recommendations")
        if hotspot_overlap > 0 or wdpa_overlap > 0:
            st.info("""
            Your farm is located in an environmentally sensitive area. Consider these sustainable practices:
            - Implement agroforestry systems to enhance biodiversity
            - Use cover crops to prevent soil erosion
            - Minimize pesticide use through integrated pest management
            - Establish wildlife corridors and habitat areas
            - Practice water conservation techniques
            """)
        else:
            st.info("""
            While your farm is not in a designated protected area or biodiversity hotspot, sustainable practices are still recommended:
            - Practice crop rotation to maintain soil health
            - Consider establishing natural borders with native plants
            - Use efficient irrigation methods
            - Monitor for invasive species
            """)
    
    # Download options
    st.sidebar.header("Download Options")
    if st.sidebar.button("Export NDVI to Google Drive"):
        task = ee.batch.Export.image.toDrive(
            image=ndvi,
            description='Farm_NDVI_Export',
            scale=10,
            region=farm_ee,
            maxPixels=1e13
        )
        task.start()
        st.sidebar.success("Export started! Check your Google Drive tasks.")

else:
    st.info("Please upload a shapefile or enter coordinates to begin the analysis.")

# Footer with info
st.markdown("---")
st.markdown("""
**Note:** This app uses Google Earth Engine data. Results are for informational purposes.  
For more details on the datasets used, visit [Google Earth Engine Data Catalog](https://developers.google.com/earth-engine/datasets/).
""")
