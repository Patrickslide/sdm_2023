import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import os
import numpy as np

world_urban = gpd.read_file("ne_10m_urban_areas")

bbox = (7, 45, 8.0, 46)

bats_bbox = np.array([[7.6, 45.4], [7.6,45.6],[7.1, 45.6],[7.1, 45.4]])

# Filter data within the bounding box
ppy_urban = world_urban.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
#ppy_cities = world_cities.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]

fig, ax = plt.subplots()

# Plot the filtered geographic data

bats_area = Polygon(bats_bbox, edgecolor="brown", facecolor="None")
xticks = list(range(45))


ppy_urban.plot(ax=ax, color='grey', edgecolor='black')
ax.set_xlim(7,7.8)
ax.add_patch(bats_area)
ax.legend(labels = ["Pygmaeus AOI"],loc="lower right")
ax.set_ylim(45.1, 45.62)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Display the map
plt.savefig('pipistrellus pygmaeus and urban areas', dpi=250, bbox_inches='tight')
plt.show()