#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[6]:


file_path = "/Users/rzhoufias.uni-frankfurt.de/Documents/PhD_Franziska/Headon/Statistical_Plots_for_Paper/Tracking data for figure 1/"

# Dataset of periderm. 
track_df_periderm_1 = pd.read_csv(file_path + "Periderm/ACTN 3.2.25 video 1 Periderm rotate 23.tif_Tracking Table Full.csv")
track_df_periderm_2 = pd.read_csv(file_path + "Periderm/ACTN 3.2.25 video 3 periderm.tif_Tracking Table Full.csv")
track_df_periderm_3 = pd.read_csv(file_path + "Periderm/ACTN 3.2.25 video 4 periderm rotate 7.tif_Tracking Table Full.csv")
track_df_periderm_4 = pd.read_csv(file_path + "Periderm/ACTN 3.2.25 video 5 Periderm rotate 7.tif_Tracking Table Full.csv")

# Dataset of epidermis. 
track_df_epidermis_1 = pd.read_csv(file_path + "Epidermis/ACTN 3.2.25 video 1 epidermis rotate 23.tif_Tracking Table Full.csv")
track_df_epidermis_2 = pd.read_csv(file_path + "Epidermis/ACTN 3.2.25 video 3 epidermis.czi_Tracking Table full.csv")
track_df_epidermis_3 = pd.read_csv(file_path + "Epidermis/ACTN 3.2.25 video 4 epidermis rotate 7.tif_Tracking Table Full.csv")
track_df_epidermis_4 = pd.read_csv(file_path + "Epidermis/ACTN 3.2.25 video 5 epidermis rotate 7.tif_Tracking Table Full.csv")

# Dataset of dermis. 
track_df_dermis_1 = pd.read_csv(file_path + "Dermis/ACTN 3.2.25 video 1 dermis rotate 23.tif_Tracking Table Full.csv")
track_df_dermis_2 = pd.read_csv(file_path + "Dermis/ACTN 3.2.25 video 3 dermis.tif_Tracking Table Full.csv")
track_df_dermis_3 = pd.read_csv(file_path + "Dermis/ACTN 3.2.25 video 4 dermis rotate 7.tif_Tracking Table Full.csv")
track_df_dermis_4 = pd.read_csv(file_path + "Dermis/ACTN 3.2.25 video 5 dermis rotate 7.tif_Tracking Table Full.csv")


# In[9]:


# Check the frame order of each data sample. 

def check_frame_order(track_df): 
    # Extract all existing track ids. 
    all_track_id = track_df["Track"].unique()
    # Get the sequence of Frames for each track id and extract the difference between two consicutive frame numbers. 
    unusual_track = [] # Store the track with unusual frames. 
    for i in all_track_id:
        frame_diff = np.diff(track_df.loc[track_df["Track"]== i, "Frame"])
        if np.all(frame_diff == 1) == False: 
            unusual_track.append(i)
    return unusual_track


# In[75]:


check_frame_order(track_df_epidermis_4)
# np.diff(track_df.loc[track_df["Track"]== track_id, "Frame"]) # Check the exact problem for single track.


# In[3]:


# Get the X and Y coordinates with the Track number.

def extract_X_Y_single_track(track_df, track_id):
    # ttrack_df: the dataframe containing track information.
    # track_num: the id of the single track. 
    track_X = track_df.loc[track_df["Track"]==track_id, "X"]
    track_Y = track_df.loc[track_df["Track"]==track_id, "Y"]
    return np.array(track_X), np.array(track_Y)


# In[4]:


# The Euclidean distance between strat and end points for each pseudo trajectory.
# Without considering the measure conversion -- in pixel. 

def euclidean_distance_pixel_single_track(x, y):
    # Calculate the Euclidean distance between the start and end points for a single track.
    start_x, end_x = x[0], x[-1]
    start_y, end_y = y[0], y[-1]
    euc_distance = np.sqrt((end_x - start_x)**2 + (end_y - start_y)**2)
    return euc_distance


# In[5]:


# Calculate the Euclidean distance in μm scale 
def euclidean_distance_single_track(track_df, track_id, pixel_conversion):
    # extract the single track x and y coordinates. 
    x, y = extract_X_Y_single_track(track_df, track_id)
    # Get the euclidean distance in pixel.
    euc_dist_pxl = euclidean_distance_pixel_single_track(x, y)
    # Convert into μm with pixel conversion factor. 
    euc_dist = euc_dist_pxl * pixel_conversion
    return euc_dist


# In[6]:


# Calculate the track length for each pseudo trajectories. 
# Without considering the measure conversion -- in pixel.

def length_pixel_single_track(x, y):
    # Calculate the differences between consecutive points for a single track.
    dx = np.diff(x)
    dy = np.diff(y)
    # Calculate the Euclidean distance for each segment and sum them
    length = np.sum(np.sqrt(dx**2 + dy**2))
    return length


# In[7]:


# Calucalte the total time for each track based on the frame numbers.

def total_frame_single_track(track_df, track_id):
    frame_num =  np.sum(np.diff(track_df.loc[track_df["Track"]==track_id, "Frame"].sort_values()))
    return frame_num


# In[8]:


# Calcualte the speed of a single track. 

def speed_single_track(track_df, track_id, time_conversion, pixel_conversion): 

    # Calculate the total time for track. 
    frame_number = total_frame_single_track(track_df, track_id)
    total_time = frame_number * time_conversion

    # Compute the total length of track.
    x, y = extract_X_Y_single_track(track_df, track_id)
    track_length = length_pixel_single_track(x, y)
    # Take account into the measure conversion. Speed = total track length / total time
    track_speed = track_length * pixel_conversion / total_time

    return track_speed


# In[9]:


# Calcualate the persistence of a single track. 
# Persistence is calcualted as the euclidean distance / total length of the track.

def persistence_single_track(track_df, track_id): 
    x, y = extract_X_Y_single_track(track_df, track_id)
    # Euclidean distance 
    euc_distance = euclidean_distance_pixel_single_track(x, y)
    # Total length of track 
    length = length_pixel_single_track(x, y)
    # Persistence calculated as euclidean distance / total track length.
    persistence = euc_distance / length
    return persistence


# In[10]:


# Add single track information to an existing data frame. 
def add_track_information(exist_dataframe, track_df, time_conversion, pixel_conversion):
    
    # Extract all existing track ids. 
    all_track_id = track_df["Track"].unique()
    # Extract the image name for each dataset table to identify cell tracks. 
    image_id = track_df["Image_ID"][0]

    for i in all_track_id:
        # Calculate the speed.
        speed = speed_single_track(track_df, i, time_conversion, pixel_conversion)
        # Calculate the persistence. 
        persistence = persistence_single_track(track_df, i)
        # Calculate the Euclidean distance. 
        euc_distance = euclidean_distance_single_track(track_df, i, pixel_conversion)
        # Contract the new row that should be add to the existing dataframe. 
        new_row = pd.DataFrame({"Image_ID": [image_id], 
                               "Track_ID": [i],
                               "Speed": [speed], 
                               "Persistence": [persistence],
                               "Euclidean_Dist": [euc_distance]})
        # Add the new row to the existing dataframe. 
        exist_dataframe = pd.concat([exist_dataframe, new_row], ignore_index=True)
        
    return exist_dataframe
    


# In[11]:


# Add tge data from a list of dataframe together in one existing dataframe. 

def summary_track_information_in_one_dataframe(df_list, time_conversion, pixel_conversion):
    final_df = pd.DataFrame(columns = ["Image_ID", "Track_ID", "Speed", "Persistence", "Euclidean_Dist"])
    for df in df_list:
        final_df  = add_track_information(final_df, df, time_conversion, pixel_conversion)
    return final_df


# In[12]:


# Dataframe for periderm. 
periderm_df_list = [track_df_periderm_1, track_df_periderm_2, track_df_periderm_3, track_df_periderm_4]
periderm_df = summary_track_information_in_one_dataframe(periderm_df_list, 10, 0.69)

# Dataframe for epidermis. 
epidermis_df_list = [track_df_epidermis_1, track_df_epidermis_2, track_df_epidermis_3, track_df_epidermis_4]
epidermis_df = summary_track_information_in_one_dataframe(epidermis_df_list, 10, 0.69)

# Dataframe for dermis. 
dermis_df_list = [track_df_dermis_1, track_df_dermis_2, track_df_dermis_3, track_df_dermis_4]
dermis_df = summary_track_information_in_one_dataframe(dermis_df_list, 10, 0.69)


# In[13]:


# Plot the speed datasets as violoin plots, grouped by cell types. 

fig = plt.figure(figsize=(2.5,4), dpi=300)

# Violin plot with log on y-scale to visualise the difference in small values better. 
plt.violinplot([periderm_df["Speed"], epidermis_df["Speed"], dermis_df["Speed"]], showmeans=True) #, showextrema=False)

# Add the mean of each single dataset. 

# For periderm data set. 
plt.scatter([1]*4, [np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 1 periderm rotate 23.tif", "Speed"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 3 periderm.tif", "Speed"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 4 periderm rotate 7.tif", "Speed"]),
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 5 periderm rotate 7.tif", "Speed"])], 
            color="orange", alpha=0.6, zorder=2)
# For epidermis data set. 
plt.scatter([2]*4, [np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 1 epidermis rotate 23.tif", "Speed"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 3 epidermis.czi", "Speed"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 4 epidermis rotate 7.tif", "Speed"]),
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 5 epidermis rotate 7.tif", "Speed"])], 
            color="orange", alpha=0.6, zorder=2)
# For dermis data set. 
plt.scatter([3]*4, [np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 1 dermis rotate 23.tif", "Speed"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 3 dermis.tif", "Speed"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 4 dermis rotate 7.tif", "Speed"]),
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 5 dermis rotate 7.tif", "Speed"])], 
            color="orange", alpha=0.6, zorder=2)

# Customize the plot. 
plt.xticks([1, 2, 3], ["Periderm", "Basal Epidermis", "Dermis"])  # Label each part
fig.autofmt_xdate(rotation=45)
#plt.yticks([-3, -2, -1, 0, 1 ])
# Hide the right and top spines
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.ylabel("Speed (μm/min)") 
# plt.savefig("speed_figure_1.svg", format="svg")


# In[1]:


# # Plot the speed datasets as violoin plots, grouped by cell types. 
# # For a better illustration of the difference in small values, applied logarithmic on y axis. 

# fig = plt.figure(figsize=(2.5,4), dpi=300)

# # Violin plot with log on y-scale to visualise the difference in small values better. 
# plt.violinplot([np.log(epidermis_df["Speed"]), np.log(periderm_df["Speed"]), np.log(dermis_df["Speed"])], showmeans=True) #, showextrema=False)

# # Add the mean of each single dataset. 
# # For epidermis data set. 
# plt.scatter([1]*3, [np.log(np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="Cell cycle chicken video 1 epidermis.tif", "Speed"])), 
#                     np.log(np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="Cell cycle chicken video 2 epidermis.tif", "Speed"])), 
#                     np.log(np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="Cell cycle chicken video 3 epidermis.tif", "Speed"]))], 
#             color="orange", alpha=0.6, zorder=2)
# # For periderm data set. 
# plt.scatter([2]*3, [np.log(np.mean(periderm_df.loc[periderm_df["Image_ID"]=="Fucci Chicken video 1 periderm projection.tif", "Speed"])), 
#                     np.log(np.mean(periderm_df.loc[periderm_df["Image_ID"]=="Fucci Chicken video 2 periderm projection.tif", "Speed"])), 
#                     np.log(np.mean(periderm_df.loc[periderm_df["Image_ID"]=="Fucci Chicken video 3 periderm projection.tif", "Speed"]))], 
#             color="orange", alpha=0.6, zorder=2)
# # For dermis data set. 
# plt.scatter([3]*3, [np.log(np.mean(dermis_df.loc[dermis_df["Image_ID"]=="Cell cycle chicken condensate 1 slice 3.tif", "Speed"])), 
#                     np.log(np.mean(dermis_df.loc[dermis_df["Image_ID"]=="Cell cycle chicken condensate 2 slice 5.tif", "Speed"])), 
#                     np.log(np.mean(dermis_df.loc[dermis_df["Image_ID"]=="Cell cycle chicken condensate 3 slice 4.tif", "Speed"]))], 
#             color="orange", alpha=0.6, zorder=2)

# # Customize the plot. 
# plt.xticks([1, 2, 3], ["Basal epidermis", "Periderm", "Dermis"])  # Label each part
# fig.autofmt_xdate(rotation=45)
# plt.yticks([-3, -2, -1, 0, 1 ])
# # Hide the right and top spines
# ax = plt.gca()
# ax.spines[['right', 'top']].set_visible(False)
# plt.ylabel("Logarithm of speed") # Speed unit (μm/min). 
# plt.savefig("log_speed_figure_1.svg", format="svg")


# In[22]:


# Plot the persistence datasets as violoin plots, grouped by cell types. 

fig = plt.figure(figsize=(2.5,4), dpi=300)

# Violin plot with log on y-scale to visualise the difference in small values better. 
plt.violinplot([periderm_df["Persistence"], epidermis_df["Persistence"], dermis_df["Persistence"]], showmeans=True) #, showextrema=False)

# Add the mean of each single dataset. 

# For periderm data set. 
plt.scatter([1]*4, [np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 1 periderm rotate 23.tif", "Persistence"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 3 periderm.tif", "Persistence"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 4 periderm rotate 7.tif", "Persistence"]),
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 5 periderm rotate 7.tif", "Persistence"])], 
            color="orange", alpha=0.6, zorder=2)
# For epidermis data set. 
plt.scatter([2]*4, [np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 1 epidermis rotate 23.tif", "Persistence"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 3 epidermis.czi", "Persistence"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 4 epidermis rotate 7.tif", "Persistence"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 5 epidermis rotate 7.tif", "Persistence"])], 
            color="orange", alpha=0.6, zorder=2)
# For dermis data set. 
plt.scatter([3]*4, [np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 1 dermis rotate 23.tif", "Persistence"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 3 dermis.tif", "Persistence"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 4 dermis rotate 7.tif", "Persistence"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 5 dermis rotate 7.tif", "Persistence"])], 
            color="orange", alpha=0.6, zorder=2)

# Customize the plot. 
plt.xticks([1, 2, 3], ["Periderm", "Basal Epidermis", "Dermis"])  # Label each part
fig.autofmt_xdate(rotation=45)
#plt.yticks([-3, -2, -1, 0, 1 ])
# Hide the right and top spines
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.ylabel("Persistence") 
# plt.savefig("persistence_figure_1.svg", format="svg")


# In[15]:


# Plot the Euclidean distance datasets as violoin plots, grouped by cell types. 

fig = plt.figure(figsize=(2.5,4), dpi=300)

# Violin plot with log on y-scale to visualise the difference in small values better. 
plt.violinplot([periderm_df["Euclidean_Dist"], epidermis_df["Euclidean_Dist"], dermis_df["Euclidean_Dist"]], showmeans=True) #, showextrema=False)

# Add the mean of each single dataset. 

# For periderm data set. 
plt.scatter([1]*4, [np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 1 periderm rotate 23.tif", "Euclidean_Dist"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 3 periderm.tif", "Euclidean_Dist"]), 
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 4 periderm rotate 7.tif", "Euclidean_Dist"]),
                    np.mean(periderm_df.loc[periderm_df["Image_ID"]=="ACTN 3.2.25 video 5 periderm rotate 7.tif", "Euclidean_Dist"])], 
            color="orange", alpha=0.6, zorder=2)
# For epidermis data set. 
plt.scatter([2]*4, [np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 1 epidermis rotate 23.tif", "Euclidean_Dist"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 3 epidermis.czi", "Euclidean_Dist"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 4 epidermis rotate 7.tif", "Euclidean_Dist"]), 
                    np.mean(epidermis_df.loc[epidermis_df["Image_ID"]=="ACTN 3.2.25 video 5 epidermis rotate 7.tif", "Euclidean_Dist"])], 
            color="orange", alpha=0.6, zorder=2)
# For dermis data set. 
plt.scatter([3]*4, [np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 1 dermis rotate 23.tif", "Euclidean_Dist"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 3 dermis.tif", "Euclidean_Dist"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 4 dermis rotate 7.tif", "Euclidean_Dist"]), 
                    np.mean(dermis_df.loc[dermis_df["Image_ID"]=="ACTN 3.2.25 video 5 dermis rotate 7.tif", "Euclidean_Dist"])], 
            color="orange", alpha=0.6, zorder=2)

# Customize the plot. 
plt.xticks([1, 2, 3], ["Periderm", "Basal Epidermis", "Dermis"])  # Label each part
fig.autofmt_xdate(rotation=45)
#plt.yticks([-3, -2, -1, 0, 1 ])
# Hide the right and top spines
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.ylabel("Euclidean Distance (μm)") 
plt.savefig("euclidean_dist_figure_1.svg", format="svg")


# In[ ]:




