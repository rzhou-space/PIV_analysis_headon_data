#!/usr/bin/env python
# coding: utf-8

# In[1]:


import h5py
import matplotlib.pyplot as plt
import numpy as np
import Polarisation
import PIV_postprocessing
from scipy.ndimage import median_filter


# # Single Aggregate Extraction 

# In[2]:


# Extract the data information.
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/subregion_cut1.h5"
with h5py.File(h5_file_path, "r") as h5f:
    green_channel = h5f["green"][:]
    red_channel = h5f["red"][:]

# Visualisation 
plt.imshow(green_channel[130])


# # Manually extract single aggregates.

# In[131]:


def extract_single_aggregate(channel):
    # Extract information from the dataset with given channel. 
    # Extract the image information for al frames. 
    aggregate = []
    for i in range(np.shape(channel)[0]):
        aggregate.append(channel[i][380:500, 765:885]) # Manually fit the region of aggregate!
    return np.array(aggregate)


# In[132]:


# Generate the subregion of images. 
aggregate = extract_single_aggregate(green_channel)
np.shape(aggregate) # (number of frames, y coordinate, x coordinate) dimensional 

# The aggregate image at different time point. 
plt.imshow(aggregate[190])


# In[20]:


# # Save the all frames of subregion into .h5 file. 
# h5_file_path = "single_aggregate_v2(2).h5"
# with h5py.File(h5_file_path, "w") as h5f:
#      h5f.create_dataset("image", data = aggregate)


# # Circle Mask 

# In[133]:


# Load data. 

# Red in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/Polarisation/single_aggregate_cut1_v2.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate = h5f["green"][:]


# In[154]:


# Function for cutting circle: circle_mask(center_x, center_y, radius, img)

circle_region = Polarisation.circle_mask(70, 50, 40, single_aggregate[150]) # x, y, and radius of circle. And the image that is cut.
plt.imshow(circle_region)


# # Calculation of the Polarisation Score

# In[2]:


# The calculation of polarisation score: 

def temperal_polarisation(img_series, center_x, center_y, radius):
    # Calculate the polarisation score for all frames (with conversion into minute) in the img_series. 
    total_frame = np.shape(img_series)[0]

    polarisation_score = np.zeros(total_frame)
    for i in range(total_frame):
        polarisation_score[i] = Polarisation.polarisation_vector(center_x[i], center_y[i], radius[i], img_series[i])[0]
    return polarisation_score


# ## Polarity and Velocity for aggregate h1.

# In[84]:


# Construction of parameter lists for calculation of polarisation score.
# Parameters for h1: 
# center_x_h1 = np.concatenate([np.repeat(46, 130), np.repeat(48, 10), np.repeat(50, 5), np.repeat(52, 5), np.repeat(56, 5), 
#                           np.repeat(60, 10), np.repeat(64, 5), np.repeat(66, 5), np.repeat(68, 5), np.repeat(70, 5), 
#                           np.repeat(72, 10), np.repeat(74, 10), np.repeat(74, 10), np.repeat(74, 5), np.repeat(74, 15), 
#                           np.repeat(74, 140)])

# center_y_h1 = np.concatenate([np.repeat(48, 130), np.repeat(50, 5), np.repeat(52, 10), np.repeat(54, 10), np.repeat(56, 5), 
#                           np.repeat(58, 5), np.repeat(56, 20), np.repeat(54, 10), np.repeat(52, 5), np.repeat(54, 5), np.repeat(56, 5), 
#                           np.repeat(60, 5), np.repeat(62, 10), np.repeat(60, 10), np.repeat(60, 140)])

center_x_h1 = np.concatenate([np.repeat(46, 120), np.repeat(48, 10), np.repeat(50, 10), np.repeat(52, 10), np.repeat(54, 10), np.repeat(56, 10), 
                             np.repeat(58, 10), np.repeat(60, 10), np.repeat(62, 10), np.repeat(64, 10), np.repeat(66, 10), 
                             np.repeat(68, 10), np.repeat(70, 125), np.repeat(68, 10), np.repeat(66, 10)])

# center_y_h1 = np.concatenate([np.repeat(46, 120), np.repeat(48, 10), np.repeat(50, 10), np.repeat(52, 10), np.repeat(54, 10), np.repeat(56, 10), 
#                              np.repeat(58, 10), np.repeat(60, 10), np.repeat(62, 10), np.repeat(64, 10), np.repeat(66, 10), 
#                              np.repeat(68, 10), np.repeat(70, 135), np.repeat(68, 10)])

radius_h1 = np.concatenate([np.repeat(40, 120), np.repeat(50, 255)]) 


# In[14]:


# Read in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/single_aggregate_h1.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate_h1 = h5f["image"][:]

# Calculate the polarisation score for the single aggregate. 
h1_polarisation = temperal_polarisation(single_aggregate_h1, center_x_h1, center_x_h1, radius_h1) / 35000 # Normalisation factor 35000


# In[86]:


# Caluclate the velocity over time for the single aggregate based on the PIV vectors. 

# Load the PIV data. 
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/piv_cut1(2)/subregion_all_frame/"
with h5py.File(h5_file_path+"U_matrix", "r") as h5f:
    U = h5f["matrix"][:] # shape: frame, rows, columns 

with h5py.File(h5_file_path+"V_matrix", "r") as h5f:
    V = h5f["matrix"][:]

with h5py.File(h5_file_path+"M_matrix", "r") as h5f:
    M = h5f["matrix"][:]

with h5py.File(h5_file_path+"xgrid_matrix", "r") as h5f:
    xgrid = h5f["matrix"][:]

with h5py.File(h5_file_path+"ygrid_matrix", "r") as h5f:
    ygrid = h5f["matrix"][:]

t0 = 0
t1 = np.shape(U)[0]
time_conversion = 15
pxl_conversion = 1.8

center_h1 = [530, 100] # center in [y, x] coordinate. 
# Average flow vector from neighbour region. 
avg_u_h1, avg_v_h1 = PIV_postprocessing.avrage_neighbour_direction_over_time(center_h1, t0, t1, 100, xgrid, ygrid, U, V)
# Get the velocity based on the average velocity vector. 
velocity_h1 = np.sqrt(avg_u_h1**2 + avg_v_h1**2) * pxl_conversion / time_conversion


# ## Polarity and Velocity for v1. 

# In[42]:


# Read in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/Polarisation/single_aggregate_cut1(2)_v1.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate_v1 = h5f["green"][:]

# Calculate the polarisation score for the single aggregate. 
v1_polarisation = Polarisation.temperal_polarisation(single_aggregate_v1, 70, 50, 50, 15) / 45000 # Normalisation factor 35000


# In[43]:


# Caluclate the velocity over time for the single aggregate based on the PIV vectors. 

# Load the PIV data. 
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/piv_cut1(2)/subregion_all_frame/"
with h5py.File(h5_file_path+"U_matrix", "r") as h5f:
    U = h5f["matrix"][:] # shape: frame, rows, columns 

with h5py.File(h5_file_path+"V_matrix", "r") as h5f:
    V = h5f["matrix"][:]

with h5py.File(h5_file_path+"M_matrix", "r") as h5f:
    M = h5f["matrix"][:]

with h5py.File(h5_file_path+"xgrid_matrix", "r") as h5f:
    xgrid = h5f["matrix"][:]

with h5py.File(h5_file_path+"ygrid_matrix", "r") as h5f:
    ygrid = h5f["matrix"][:]

t0 = 0
t1 = np.shape(U)[0]
time_conversion = 15
pxl_conversion = 1.8

center_v1 = [450, 690] # center in [y, x] coordinate. 
# Average flow vector from neighbour region. 
avg_u_v1, avg_v_v1 = PIV_postprocessing.avrage_neighbour_direction_over_time(center_v1, t0, t1, 100, xgrid, ygrid, U, V)
# Get the velocity based on the average velocity vector. 
velocity_v1 = np.sqrt(avg_u_v1**2 + avg_v_v1**2) * pxl_conversion / time_conversion


# ## Polarity and Velocity for h2. 

# In[10]:


# center_x_h2 = np.concatenate([np.repeat(20, 135), np.repeat(22, 5), np.repeat(24, 5), np.repeat(26, 5), np.repeat(28, 5), 
#                              np.repeat(30, 5), np.repeat(30, 5), np.repeat(32, 5), np.repeat(36, 5), np.repeat(40, 5), 
#                              np.repeat(44, 5), np.repeat(46, 5), np.repeat(50, 5), np.repeat(52, 5), np.repeat(56, 5), 
#                              np.repeat(60, 5), np.repeat(62, 5), np.repeat(58, 5), np.repeat(60, 160)])

center_x_h2 = np.concatenate([np.repeat(20+25, 135), np.repeat(22+25, 5), np.repeat(24+25, 5), np.repeat(26+25, 5), np.repeat(28+25, 5), 
                             np.repeat(30+25, 5), np.repeat(30+25, 5), np.repeat(32+25, 5), np.repeat(36+25, 5), np.repeat(40+25, 5), 
                             np.repeat(44+25, 5), np.repeat(46+25, 5), np.repeat(50+25, 5), np.repeat(52+15, 5), np.repeat(56+15, 5), 
                             np.repeat(60+15, 5), np.repeat(62+15, 5), np.repeat(58+15, 5), np.repeat(60+15, 160)])

center_y_h2 = np.concatenate([np.repeat(48, 135), np.repeat(46, 5), np.repeat(46, 5), np.repeat(46, 5), np.repeat(44, 5), 
                             np.repeat(46, 5), np.repeat(50, 5), np.repeat(52, 5), np.repeat(52, 5), np.repeat(52, 5), 
                             np.repeat(52, 5), np.repeat(54, 5), np.repeat(54, 5), np.repeat(56, 5), np.repeat(56, 5), 
                             np.repeat(56, 5), np.repeat(54, 5), np.repeat(54, 5), np.repeat(56, 160)])

radius_h2 = np.concatenate([np.repeat(40, 180), np.repeat(40, 20), np.repeat(40, 175)])


# In[11]:


# Read in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "single_aggregate_h2(2).h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate_h2 = h5f["image"][:]

# Calculate the polarisation score for the single aggregate. 
h2_polarisation = temperal_polarisation(single_aggregate_h2, center_x_h2, center_x_h2, radius_h2) / 15000 # Normalisation factor 35000


# In[12]:


# Caluclate the velocity over time for the single aggregate based on the PIV vectors. 

# Load the PIV data. 
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/piv_cut1/frame_1_to_375/"
with h5py.File(h5_file_path+"U_matrix", "r") as h5f:
    U = h5f["matrix"][:] # shape: frame, rows, columns 

with h5py.File(h5_file_path+"V_matrix", "r") as h5f:
    V = h5f["matrix"][:]

with h5py.File(h5_file_path+"M_matrix", "r") as h5f:
    M = h5f["matrix"][:]

with h5py.File(h5_file_path+"xgrid_matrix", "r") as h5f:
    xgrid = h5f["matrix"][:]

with h5py.File(h5_file_path+"ygrid_matrix", "r") as h5f:
    ygrid = h5f["matrix"][:]

t0 = 0
t1 = np.shape(U)[0]
time_conversion = 15
pxl_conversion = 1.8

center_h2 = [190, 170] # center in [y, x] coordinate. 
# Average flow vector from neighbour region. 
avg_u_h2, avg_v_h2 = PIV_postprocessing.avrage_neighbour_direction_over_time(center_h2, t0, t1, 100, xgrid, ygrid, U, V)
# Get the velocity based on the average velocity vector. 
velocity_h2 = np.sqrt(avg_u_h2**2 + avg_v_h2**2) * pxl_conversion / time_conversion


# ## Polarity and Velocity for v2.

# In[35]:


# center_x_v2 = np.concatenate([np.repeat(95, 120), np.repeat(92, 10), np.repeat(90, 10), np.repeat(88, 10), np.repeat(86, 10),
#                               np.repeat(84, 10), np.repeat(82, 10), np.repeat(80, 10), np.repeat(78, 10), np.repeat(76, 10), 
#                               np.repeat(74, 10), np.repeat(72, 10), np.repeat(70, 10), np.repeat(68, 10), np.repeat(66, 10),
#                               np.repeat(64, 10), np.repeat(62, 10), np.repeat(60, 95)])
# center_y_v2 = np.concatenate([np.repeat(45, 190), np.repeat(50, 185)]) 
# radius_v2 = np.concat([np.repeat(40, 190), np.repeat(50, 185)])

# center_x_v2 = np.repeat(60, 375)
# center_y_v2 = np.repeat(50, 375)
# radius_v2 = np.repeat(50, 375)


# In[36]:


# Read in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/Polarisation/single_aggregate_cut1_v2.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate_v2 = h5f["green"][:]

# Calculate the polarisation score for the single aggregate. 
# v2_polarisation = temperal_polarisation(single_aggregate_v2, center_x_v2, center_y_v2, radius_v2) / 30000 # Normalisation factor 30000
v2_polarisation = Polarisation.temperal_polarisation(single_aggregate_v2, 60, 50, 50, 15) / 30000


# In[25]:


# Caluclate the velocity over time for the single aggregate based on the PIV vectors. 

# Load the PIV data. 
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/piv_cut1/frame_1_to_375/"
with h5py.File(h5_file_path+"U_matrix", "r") as h5f:
    U = h5f["matrix"][:] # shape: frame, rows, columns 

with h5py.File(h5_file_path+"V_matrix", "r") as h5f:
    V = h5f["matrix"][:]

with h5py.File(h5_file_path+"M_matrix", "r") as h5f:
    M = h5f["matrix"][:]

with h5py.File(h5_file_path+"xgrid_matrix", "r") as h5f:
    xgrid = h5f["matrix"][:]

with h5py.File(h5_file_path+"ygrid_matrix", "r") as h5f:
    ygrid = h5f["matrix"][:]

t0 = 0
t1 = np.shape(U)[0]
time_conversion = 15
pxl_conversion = 1.8

center_v2 = [440, 825] # center in [y, x] coordinate. 
# Average flow vector from neighbour region. 
avg_u_v2, avg_v_v2 = PIV_postprocessing.avrage_neighbour_direction_over_time(center_v2, t0, t1, 100, xgrid, ygrid, U, V)
# Get the velocity based on the average velocity vector. 
velocity_v2 = np.sqrt(avg_u_v2**2 + avg_v_v2**2) * pxl_conversion / time_conversion


# # Plot the Twinplot for velocity and polarisation. 

# In[21]:


# Plot the temporal polarisation score and the velocity together over time 
# polarisation score for all frames (with conversion into minute) in the img_series. 

def twinplot_polarisation_velocity(avg_u, time_conversion, velocity, polarisation, img_name):
    # Plot the normalised polarisation score.  
    fig, ax = plt.subplots(1,1, facecolor=(1, 1, 1), dpi=300, figsize=(4, 2))
    # Plot with two y-axis.
    twin_stacked = ax.twinx()
    # Convert frame into minutes with time_conversion.
    time_axis = np.array([(i * time_conversion) for i in range(len(avg_u))])
    
    # Plot the line for the velocity socre. 
    # Plot the results in the background. 
    p1 = ax.scatter(time_axis, velocity, alpha=0.2, c="orange", zorder=0)
    # Plot the median filter at the top. 
    p2 = ax.plot(time_axis, median_filter(velocity, size=30, mode="nearest"), c="orange", zorder=1, label="Velocity")
    ax.set_ylabel('Velocity (μm/min)')
    
    # Plot for polarisation. 
    twin1 = twin_stacked.plot(time_axis, polarisation[5:370], alpha=0.3, zorder=2, c="blue")
    # Apply median filter for smoothing.
    twin2 = twin_stacked.plot(time_axis, median_filter(polarisation[5:370], size=50, mode="nearest"), c="blue", zorder=3, label="Polarisation")
    twin_stacked.set_ylabel("Polarisation")
    
    # Add line to highlight the maximum velocity (based on median filtered).
    maxi_velocity_frame = np.argmax(median_filter(velocity, size=30, mode="nearest"))
    ax.axvline(x = maxi_velocity_frame * time_conversion, color="red", zorder=4)
    
    ax.set_xlabel("Time (min) ")
    plt.savefig(img_name, format="svg")
    fig.legend(bbox_to_anchor=(.95, 1.17))
    plt.show()


# In[45]:


twinplot_polarisation_velocity(avg_u_v1, 15, velocity_v1, v1_polarisation, "v1_polarisation_velocity.svg")


# In[ ]:




