#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author Vladimir S. FONOV & Sofia FERNANDEZ-LOZANO
# @date 13/11/2024

import argparse
import warnings
import re
import torch
from minc2_simple import minc2_file
import torch.nn.functional as F
from torch.serialization import SourceChangeWarning
#TODO: change cnn_subrepo name
from model.util import segment_with_patches_overlap


def load_minc_volume(file_path):
    """
    Load a MINC file as a PyTorch tensor.
    """
    minc_vol = minc2_file(file_path)
    minc_vol.setup_standard_order()
    return minc_vol.load_complete_volume_tensor(minc2_file.MINC2_FLOAT)


def crop_or_pad_volume(volume, params):
    """
    Apply cropping or padding to the volume based on user-specified parameters.
    """
    if params.cropvol > 0:
        volume = volume[
            :,
            :,
            params.cropvol:-params.cropvol,
            params.cropvol:-params.cropvol,
            params.cropvol:-params.cropvol
        ]
    elif params.padvol > 0:
        volume = F.pad(
            volume,
            (params.padvol,) * 6,
            mode='constant',
            value=params.padfill
        )
    return volume


def segment_volume_with_model(volume, model, params):
    """
    Segment the volume using the pre-trained model with optional fuzzy output.
    """
    return segment_with_patches_overlap(
        volume,
        model,
        patch_sz=params.patch_sz,
        use_cuda=not params.cpu,
        crop=params.crop,
        bck=params.bck,
        stride=params.stride,
        out_fuzzy=params.fuzzy
    )


def apply_mask(volume, mask_path, background_label):
    """
    Apply a mask to the segmented volume, setting regions outside the mask
    to the background label.
    """
    mask = load_minc_volume(mask_path)
    volume[mask < 1] = background_label


def save_minc_volume(volume, output_path, params):
    """
    Save a PyTorch tensor volume to a MINC file.
    """
    minc_vol = minc2_file(params.in_img)  # Using the input file for metadata
    out_file = minc2_file()
    out_file.define(
        minc_vol.store_dims(),
        minc2_file.MINC2_BYTE,
        minc2_file.MINC2_BYTE
    )
    out_file.create(output_path)
    out_file.setup_standard_order()
    out_file.copy_metadata(minc_vol)
    out_file.save_complete_volume_tensor(volume.byte().contiguous())


def apply_cnn_segmentation(params):
    """
    Apply the CNN model to segment the input minc volume
    and save the result to an output file.
    """
    # Load the pre-trained model
    model = torch.load(params.model)
    if params.cpu:
        model.cpu()
    else:
        model.cuda()
    model.eval()

    # Parse input specification
    input_data = []
    input_spec_match = re.match(r"\[(.*)\]", params.in_img)
    if input_spec_match:
        # Handle multiple inputs defined by constants or files
        for element in input_spec_match[1].split(","):
            if re.match(r"^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$", element):
                # If it's a number
                input_data.append(float(element))
            else:
                # If it's a file path
                file_data = load_minc_volume(element)
                input_data.append(file_data)
        # Stack inputs along the channel dimension
        input_tensor = torch.cat(
            [torch.full(input_data[0].shape, v)
             if isinstance(v, float)
             else v for v in input_data],
            dim=1
        )
    else:
        # Handle single input file
        input_tensor = load_minc_volume(params.in_img).unsqueeze(0).unsqueeze(0)

        # Additional channels
        if params.add:
            channels = [input_tensor] + \
                [load_minc_volume(f).unsqueeze(0).unsqueeze(0)
                 for f in params.add]
            input_tensor = torch.cat(channels, dim=1)

        # Dummy channels if specified
        if params.channels > 1:
            dummy_channel = torch.full_like(input_tensor, params.fill)
            input_tensor = torch.cat(
                [input_tensor] + \
                [dummy_channel] * \
                (params.channels - 1),
                dim=1
            )

    # Handle cropping and padding
    input_tensor = crop_or_pad_volume(input_tensor, params)

    # Run segmentation
    if params.fuzzy:
        segment_result, fuzzy_result = segment_volume_with_model(
            input_tensor, model, params
        )
    else:
        segment_result = segment_volume_with_model(
            input_tensor, model, params
        )

    # Post-processing: masking, padding, and saving the result
    if params.mask:
        apply_mask(segment_result, params.mask, params.bck)

    save_minc_volume(segment_result, params.output, params)


# Main script execution
if __name__ == '__main__':
    def parse_options():
        """
        Parse command-line options for applying a pre-trained
        CNN model to segment a volume.
        """
        parser = argparse.ArgumentParser(
            description='Apply a pre-trained model to segment a MINC file',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        # Model and file inputs
        parser.add_argument(
            "model",
            type=str,
            help="Path to the pre-trained model file"
        )

        parser.add_argument(
            "in_img",
            type=str,
            help="Input MINC file or input spec in [a,b,...] format"
        )

        parser.add_argument(
            "output",
            type=str,
            help="Output MINC file"
        )

        # Optional parameters
        parser.add_argument(
            "--add",
            type=str,
            nargs='+',
            help="Additional input MINC files for extra channels"
        )
        parser.add_argument(
            "--patch_sz",
            type=int,
            default=64,
            help="Patch size for segmentation"
        )
        parser.add_argument(
            "--stride",
            type=int,
            default=None,
            help="Stride for overlapping patches (default: patch_sz-crop*2)"
        )
        parser.add_argument(
            "--channels",
            type=int,
            default=1,
            help="Number of input channels, fills with constant if >1"
        )
        parser.add_argument(
            "--crop",
            type=int,
            default=0,
            help="Crop edges of patch for overlapping segmentation"
        )
        parser.add_argument(
            "--cropvol",
            type=int,
            default=0,
            help="Crop edges of the input volume before segmentation"
        )
        parser.add_argument(
            "--padvol",
            type=int,
            default=0,
            help="Pad the input volume before segmentation"
        )
        parser.add_argument(
            "--padfill",
            type=float,
            default=0, help="Padding value for extended areas"
        )
        parser.add_argument(
            "--mask",
            type=str,
            help="Mask file to apply on the result"
        )
        parser.add_argument(
            "--bck",
            type=int,
            default=0,
            help="Background label for masked areas"
        )
        parser.add_argument(
            "--cpu",
            action="store_true",
            help="Run segmentation on CPU"
        )
        parser.add_argument(
            "-q", "--quiet",
            action="store_true",
            help="Suppress warnings"
        )
        parser.add_argument(
            "-F", "--fuzzy",
            action="store_true",
            help="Generate fuzzy output volume(s)"
        )
        parser.add_argument(
            "--fill",
            type=float,
            default=38.81240207,
            help="Fill value for missing channels"
        )

        return parser.parse_args()

    # Parse command-line arguments and suppress warnings if required
    PARAMS = parse_options()
    if PARAMS.quiet:
        warnings.filterwarnings("ignore", category=SourceChangeWarning)

    # Run segmentation
    apply_cnn_segmentation(params=PARAMS)
