#!/usr/bin/env python
# coding: utf-8

# In[1]:


import h5py
import matplotlib.pyplot as plt
import numpy as np
import Polarisation
from scipy.ndimage import median_filter


# # Single Aggregate Extraction 

# In[2]:


# Extract the data information.
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/subregion_cut1(2).h5"
with h5py.File(h5_file_path, "r") as h5f:
    green_channel = h5f["green"][:]
    red_channel = h5f["red"][:]

# Visualisation 
plt.imshow(green_channel[130])


# # Manually extract single aggregates.

# In[3]:


def extract_single_aggregate(channel):
    # Extract information from the dataset with given channel. 
    # Extract the image information for al frames. 
    aggregate = []
    for i in range(np.shape(channel)[0]):
        aggregate.append(channel[i][460:590, 40:160]) # Manually fit the region of aggregate!
    return np.array(aggregate)


# In[4]:


# Generate the subregion of images. 
aggregate = extract_single_aggregate(green_channel)
np.shape(aggregate) # (number of frames, y coordinate, x coordinate) dimensional 

# The aggregate image at different time point. 
plt.imshow(aggregate[150])


# In[5]:


# Save the all frames of subregion into .h5 file. 
h5_file_path = "single_aggregate_h2.h5"
with h5py.File(h5_file_path, "w") as h5f:
     h5f.create_dataset("image", data = aggregate)


# # Circle Mask 

# In[6]:


# Load data. 

# Red in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/single_aggregate_h1.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate = h5f["image"][:]


# In[7]:


def circle_mask(center_x, center_y, radius, img):
    # Extract a circle on the image with given center and radius.
    # center_x, center_y: integers. Define the center pixel of circle.
    # radius: integer. 
    
    # Create a grid of coordinates
    height, width = np.shape(img)[:2]
    Y, X = np.ogrid[:height, :width]
    
    # Calculate the mask for the circle
    distance_from_center = np.sqrt((X - center_x)**2 + (Y - center_y)**2)
    circular_mask = distance_from_center <= radius
    
    # Apply the mask to get pixels within the circle
    circular_region = np.zeros_like(img)
    circular_region[circular_mask] = img[circular_mask]
    
    return circular_region


# In[57]:


# Function for cutting circle: circle_mask(center_x, center_y, radius, img)

circle_region = circle_mask(48, 48, 40, single_aggregate[130]) # x, y, and radius of circle. And the image that is cut.
plt.imshow(circle_region)


# # Calculation of the Polarisation Score

# In[10]:


# The calculation of polarisation score: 

def temperal_polarisation(img_series, center_x, center_y, radius):
    # Calculate the polarisation score for all frames (with conversion into minute) in the img_series. 
    total_frame = np.shape(img_series)[0]

    polarisation_score = np.zeros(total_frame)
    for i in range(total_frame):
        polarisation_score[i] = Polarisation.polarisation_vector(center_x[i], center_y[i], radius[i], img_series[i])[0]
    return polarisation_score


# In[61]:


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


# In[62]:


# Red in the single aggregate subregion image. (could for example be generated with capture "single aggregate extraction").
h5_file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/RedBeads/single_aggregate_h1.h5"
with h5py.File(h5_file_path, "r") as h5f:
    single_aggregate_h1 = h5f["image"][:]


h1_polarisation = temperal_polarisation(single_aggregate_h1, center_x_h1, center_x_h1, radius_h1) / 35000 # Normalisation factor 35000


# In[63]:


# Plot the temporal polarisation score over time 
# polarisation score for all frames (with conversion into minute) in the img_series. 
total_frame = np.shape(single_aggregate_h1)[0]

# Plot the normalised polarisation score.  
plt.figure(figsize=(5, 2), dpi=300)

# Convert frame into minutes with time_conversion.
frame_to_min = np.array([i for i in range(total_frame)]) #*time_conversion

# Plot the results in the background. 
plt.plot(frame_to_min, h1_polarisation, c = "blue", alpha=0.3, zorder=0)
# Plot the median filter at the top. 
plt.plot(frame_to_min, median_filter(h1_polarisation, size=50, mode="nearest"), c="blue", zorder=1)

# Add verticle line: highlight the frame where the velocity is maximal (from RedBeads/PIV_postprocessing)
plt.axvline(x = 140, color="red")

plt.xlabel("frame")
plt.ylabel("polarisation")
plt.savefig("h1_polarisation_score_v1.svg", format="svg")
plt.show()


# In[ ]:




