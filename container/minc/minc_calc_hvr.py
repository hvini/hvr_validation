#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author Sofia Fernandez-Lozano
# @date 03/12/2024
import argparse
from pathlib import Path
from warnings import warn
from minc.minc_tools import mincTools


def extract_label_volumes(minc_path):
    """
    Use mincTools to print all volumes of labels contained in a MINC volume.

    Args:
        minc_path (Path): Path to the segmented image file.

    Returns:
        Dictionary: Volume in voxel-size for each label.

    Raises:
        FileNotFoundError: If the file does not exist.
    """

    if not minc_path.exists():
        raise FileNotFoundError(
            f"MINC volume '{str(minc_path)}' does not exist."
        )

    # Run print_all_labels to extract labels' volumes
    with mincTools() as minc:
        output = minc.execute_w_output([
            "print_all_labels",
            str(minc_path)
        ])

    vols_dict = {}

    for line in output.strip().split('\n'):
        # Split each line by spaces to isolate label and value
        parts = line.split()
        if len(parts) == 3 and parts[0] == 'Label:':
            label = int(parts[1])
            value = int(parts[2])
            vols_dict[label] = value

    return vols_dict


def map_labels_to_rois(label_volumes, label_mapping):
    """
    Recursively map label volumes to HC and VC for left/right hemispheres.

    Args:
        label_volumes (dict):
            Dictionary of {label: volume}.
        label_mapping (dict):
            Nested dictionary specifying which labels correspond to
            Left/Right HC/VC.
            if bilateral dict is nested — {"Left/Right": {"HC/VC": label}},
            otherwise - {"HC/VC": label}

    Returns:
        dict:
            A dictionary waith the same structure as label_mapping containing
            volume data.

    Raises:
        TypeError: If label_volumes or label_mapping are not dictionaries.
    """
    def recursive_mapping(mapping):
        """
        Helper function for recursive mapping.
        Args: mapping (dict): Current level of the label_mapping.
        Returns: dict: A nested dictionary with volumes at the leaves.
        """
        result = {}
        for key, value in mapping.items():
            if isinstance(value, dict):
                # Recur if the value is a nested dictionary
                result[key] = recursive_mapping(value)
            else:
                # Otherwise, replace the label with its corresponding volume
                result[key] = label_volumes.get(value, 0)
        return result

    if not isinstance(label_volumes, dict):
        raise TypeError("Argument 'label_volumes' is not a dictionary.")

    if not isinstance(label_mapping, dict):
        raise TypeError("Argument 'label_mapping' is not a dictionary.")

    return recursive_mapping(label_mapping)


def calculate_hvr(hc_vol, vc_vol):
    """
    Calculate HVR.

    Args:
        hc_vol (int): Volume for hippocampus.
        vc_vol (int): Volume for temporal horn of lateral ventricle.

    Returns:
        float: Hippocampal-to-Ventricle ratio.

    Raises:
        TypeError: If hc_vol or vc_vol are not integers or floats.
        ValueError: If both hc_vol and vc_vol are 0 (to avoid division by 0).
    """
    # Sanity checks
    if not isinstance(hc_vol, (int, float)):
        raise TypeError(
            f"hc_vol must be an int or float, got {type(hc_volume).__name__}"
        )

    if not isinstance(vc_vol, (int, float)):
        raise TypeError(
            f"vc_vol must be an int or float, got {type(vc_vol).__name__}"
        )

    # Avoid division by zero
    total_vol = hc_vol + vc_vol
    if total_vol == 0:
        raise ValueError(
            "The sum of HC and VC volumes cannot be zero (division by zero)."
        )

    # Calculate HVR
    hvr = hc_vol / total_vol
    return hvr


def extract_hvr_from_img(
    minc_path,
    label_mapping,
    bilateral=False,
    hc_key="HC",
    vc_key="VC",
    left_key="Left",
    right_key="Right",
):
    """
    Wrapper function to calculate HVR, by hemispheres or not, from a segmented
    MINC labels volume.

    Args:
        minc_path (str or Path): Path to the MINC volume.
        label_mapping (dict): Mapping of labels to HC/VC ROIs.
        bilateral (bool): Whether label_mapping contains laterality data.
        hc_key: Key for the Hippocampus.
        vc_key: Key for the Ventricle.
        left_key: Key for the left hemisphere.
        right_key: Key for the right hemisphere.

    Returns:
        dict: HVR values —
            {HVR: value} if !bilateral; {HVR: {L/R: value}} if bilateral.
    """
    # Extract label volumes
    minc_path = Path(minc_path)
    volumes = extract_label_volumes(minc_path)

    # Map labels to ROIs
    roi_volumes = map_labels_to_rois(volumes, label_mapping)

    # Calculate HVR
    def get_hvr_by_hemisphere(hc_key, vc_key, volume_dict):
        """
        Helper function to extract HVR by a single hemisphere.

        Args:
            hc_key (str): Key for Hippocampus.
            vc_key (str): Key for Ventricle.
            volume_dict (dict): Dictionary with the volume data.

        Returns:
            float: HVR value.
        """
        for roi, key in {"Hippocampus": hc_key, "Ventricle": vc_key}.items():
            if key not in volume_dict:
                raise KeyError(
                    f"ROI key {key} for {roi} was not found in the dictionary "
                    "for label_mapping."
                )
        hvr = calculate_hvr(
            hc_vol=volume_dict[hc_key],
            vc_vol=volume_dict[vc_key]
        )

        return hvr

    output = {}
    if bilateral:
        for side, key in {"Left": left_key, "Right": right_key}.items():
            if key not in roi_volumes:
                raise KeyError(
                    "Argument bilateral was set to True, but Key "
                    f"{key} for {side} hemisphere was not found in "
                    "the dictionary for label_mapping."
                )
            hvr = get_hvr_by_hemisphere(hc_key, vc_key, roi_volumes[key])
            output[side] = {
                "HC": roi_volumes[key].get(hc_key),
                "VC": roi_volumes[key].get(vc_key),
                "HVR": hvr,
            }
    else:
        output[side] = {
            "HC": roi_volumes.get(hc_key),
            "VC": roi_volumes.get(vc_key),
            "HVR": get_hvr_by_hemisphere(hc_key, vc_key, roi_volumes)
        }

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate HVR from a segmented MINC volume."
    )
    parser.add_argument("minc_volume", help="Path to the segmented MINC file.")
    parser.add_argument(
        "--labels",
        nargs="+",
        type=int,
        required=True,
        metavar=("HC", "VC", "[Right_HC]", "[Right_VC]"),
        help=(
            "Labels for the regions of interest. "
            "Provide 2 labels for HC and VC, or 4 labels "
            "for Left HC, Left VC, Right HC, and Right VC."
        )
    )
    args = parser.parse_args()

    if len(args.labels) == 2:
        bilateral = False
        label_mapping = {"HC": args.labels[0], "VC": args.labels[1]}
    elif len(args.labels) == 4:
        bilateral = True
        label_mapping = {
            "Left": {"HC": args.labels[0], "VC": args.labels[1]},
            "Right": {"HC": args.labels[2], "VC": args.labels[3]},
        }
    else:
        raise ValueError(
            "Either 2 (HC and VC) or 4 (L_HC, L_VC, R_HC, R_VC) labels "
            "must be provided."
        )

    hvr_results = extract_hvr_from_img(
        Path(args.minc_volume),
        label_mapping,
        bilateral,
    )

    # Print results
    print(f"Hippocampal-to-Ventricle Ratio results for {args.minc_volume}:")
    if bilateral:
        print(f"    Left hemisphere HVR: {hvr_results['HVR']['Left']:.3f}")
        print(f"    Right hemisphere HVR: {hvr_results['HVR']['Right']:.3f}")
    else:
        print(f"    HVR: {hvr_results['HVR']:.3f}")
