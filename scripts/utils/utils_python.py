#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:54:00 2024

@author: mlurig@ad.ufl.edu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:52:50 2024

@author: mlurig@ad.ufl.edu
"""
#%% imports
import copy
import cv2
from contextlib import redirect_stdout
from io import StringIO
import numpy as np
import torch
from torchvision import transforms
import phenopype as pp
import unicom

from phenopype.core.segmentation import detect_contour


#%% functions


class WarpModule(torch.nn.Module):
    def __init__(self, model) -> None:
        super().__init__()
        self.model = model

    def forward(self, x):
        return self.model(x)

def setup_unicom(model_name, resize=True):
    # Load the model from a specified source
    model, _ = unicom.model.load(model_name)
    model.cuda()
    model.eval()

    # Set up preprocessing
    pipeline = [
        transforms.ToTensor(),
        transforms.Normalize((0.48145466, 0.4578275, 0.40821073),
                             (0.26862954, 0.26130258, 0.27577711)),
        ]  
  
    ## preprend resizing step if necessary
    if resize:
        pipeline = [
            transforms.Resize(336, interpolation=transforms.functional.InterpolationMode.BICUBIC),
            transforms.CenterCrop(336),
            ] + pipeline
    
  
    ## final pipeline
    preprocessing = transforms.Compose(pipeline)

    # Wrap the model using WarpModule
    wrapped_model = WarpModule(model)

    return wrapped_model, preprocessing


def filter_mask(image, mask, min_area, ret_area=False):    

        ## convert to binary mask
        mask = mask.astype(np.uint8) * 255

        ## detect all contours above threshold
        with redirect_stdout(StringIO()):
            annotations = detect_contour(mask, min_area=min_area, stats_mode="circle")
            
        if annotations["contour"]["a"]["data"]["n"] > 0:
            countour_coords = annotations["contour"]["a"]["data"]["contour"]
            contour_info = annotations["contour"]["a"]["data"]["support"]

            ## find largest contour
            idx_largest = max(range(len(contour_info)), key=lambda i: contour_info[i]['area']) 
            coords = countour_coords[idx_largest]
            contour_info = contour_info[idx_largest]
            contour_info = {key: value for key, value in contour_info.items() if not key[:9] == "hierarchy"}

            ## draw mask
            rx, ry, rw, rh = cv2.boundingRect(coords)
            object_mask = np.zeros((rh, rw), dtype="uint8")
            object_mask = cv2.drawContours(
                image=object_mask,
                contours=[coords],
                contourIdx=0,
                thickness=-1,
                color=255,
                offset=(-rx, -ry),
            )          
            
            ## extract image and change background to white 
            object_image = copy.deepcopy(image[ry : ry + rh, rx : rx + rw])
            object_image[object_mask==0] = 0
            
            # Convert to RGBA
            rgba_image = cv2.cvtColor(object_image, cv2.COLOR_RGB2RGBA)
            rgba_image[:, :, 3] = object_mask
            
            ## prep info
            info = {
                "area": contour_info["area"],
                "bbox":  [rx, ry, rw, rh],
                "center": contour_info["center"],
                "diameter": contour_info["diameter"],
                }
            
            return rgba_image, info
        
        else:
            return None, {}

def channel_parsimony(pixels, mask, cutoffs=[0.25, 0.5, 0.75, 0.95]):
    """
    Calculates the number of unique pixel color assignments needed to reach given cutoff proportions
    of total pixel counts. Additionally, calculates mean, variance, skewness, and CV (coefficient of variation)
    on the sorted pixel counts.

    Args:
        pixels (numpy.ndarray): Array of pixel values.
        mask (numpy.ndarray): Array indicating which pixels are relevant (non-zero values in the mask).
        cutoffs (list): List of proportions of total pixel counts to cover. Default is [0.25, 0.5, 0.75, 0.95].

    Returns:
        dict: A flat dictionary containing:
            - 'mean': Mean of sorted pixel counts.
            - 'variance': Variance of sorted pixel counts.
            - 'bin_cutoffX': Number of bins needed for each cutoff (where X is the cutoff value).
    """
    # Get unique pixel values and their counts
    bins, counts = np.unique(pixels[mask > 0], return_counts=True)  # Apply mask

    # Sort counts in descending order
    sorted_indices = np.argsort(-counts)
    sorted_counts = counts[sorted_indices]
    
    # Calculate total pixel count
    total_pixels = np.sum(sorted_counts)
    
    # Calculate proportional counts and cumulative sums
    prop_counts = sorted_counts / total_pixels
    cumulative_props = np.cumsum(prop_counts)

    # Initialize the results dictionary
    results = {
        'binvar': np.var(prop_counts),
    }

    # For each cutoff, find the number of bins needed and add them to the flat dictionary
    for cutoff in cutoffs:
        num_bins_needed = np.searchsorted(cumulative_props, cutoff)
        results[f'bin{cutoff}'] = num_bins_needed

    return results

def channel_to_bins(channel, mask=None, max_val=255, blur=5, n_bins=10):
    """
    Splits a channel into n bins, applies blurring, and masks the non-relevant areas.
    Each bin is replaced by the median of the values that fall into that bin.

    Args:
        channel (numpy.ndarray): 2D array representing the channel (e.g., hue, saturation).
        mask (numpy.ndarray): Binary mask indicating relevant areas (0 for ignore). Optional.
        max_val (int): The maximum value of the channel (default 255 for 8-bit channels, 180 for hue).
        blur (int): Strength of the blur to apply. Default is 9.
        n_bins (int): Number of bins to split the channel into. Default is 18.

    Returns:
        numpy.ndarray: 2D array where each value is replaced by the median of its bin.
    """
    # Define the bin edges for splitting the channel
    bin_edges = np.linspace(0, max_val, n_bins + 1)

    # Apply blurring to the channel
    image_blurred = pp.preprocessing.blur(channel, blur)

    # Apply mask to ignore irrelevant areas if provided
    if mask is not None:
        image_blurred[mask == 0] = 0  # Set irrelevant areas to 0

    # Assign each pixel to a bin
    binned_channel = np.digitize(image_blurred, bin_edges) - 1  # Subtract 1 for 0-based indexing

    # Calculate the median for each bin
    output_channel = np.copy(image_blurred)  # Initialize output
    for i in range(n_bins):
        # Find the indices of the pixels that belong to the current bin
        bin_mask = binned_channel == i
        if np.any(bin_mask):
            # Calculate the median of the pixel values in this bin
            bin_median = np.median(image_blurred[bin_mask])
            # Assign the median to all pixels in this bin
            output_channel[bin_mask] = bin_median

    return output_channel.astype(np.uint8)