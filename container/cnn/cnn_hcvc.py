#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author Sofia Fernandez-Lozano
# @date 14/11/2024

from pathlib import Path
from argparse import Namespace
from warnings import warn
from pprint import pformat
from minc.minc_tools import mincTools
from minc.minc_mtl_qc import MincQC
from cnn.apply_segm_model import apply_cnn_segmentation
from app import PROJECT_DIR


class CNN_HCVC:
    def __init__(
        self,
        model,
        in_img=None,
        output=None,
        sep_by_sides=True,
        reference_img=None,
        work_dir=None,
        clobber=False,
        remap=None,
        add=None,
        patch_sz=96,
        stride=32,
        channels=1,
        crop=8,
        cropvol=0,
        padvol=0,
        padfill=0,
        mask=None,
        bck=0,
        cpu=True,
        quiet=False,
        fuzzy=False,
        fill=38.81240207,
    ):
        """
        Initialize the HVR CNN class.

        Args:
            model (str or Path):
                Predefined model or path to the trained CNN model.
            in_img (str or Path, optional):
                Path to the input image file. Can be updated later.
            output (str or Path, optional):
                Path to save the output file. Can be updated later.
            sep_by_sides (bool):
                Whether the input image needs to be run by hemispheres.
            reference_img (str or Path, optional if model is predefined):
                Path to the reference image used for the model application.
            work_dir (str or Path, optional):
                Path to the working directory.
            kwargs: Additional parameters to override or extend.
        """
        self.weights, self.ref_img, self.remap = self._validate_model(
            model, reference_img, remap
        )

        self.model = model.lower() \
            if model.lower() in ("simple", "detailed") else "custom"
        self.in_img = Path(in_img) if in_img else None
        self.output = Path(output) if output else None
        self.work_dir = Path(work_dir) if work_dir else None
        self.clobber = clobber
        self.sep_by_sides = sep_by_sides
        self.has_run = False

        # Initialize Params object
        self.params = Namespace(
            model=str(self.weights),
            in_img=str(self.in_img) if self.in_img else None,
            output=str(self.output) if self.output else None,
            add=add,
            patch_sz=patch_sz,
            stride=stride,
            channels=channels,
            crop=crop,
            cropvol=cropvol,
            padvol=padvol,
            padfill=padfill,
            mask=mask,
            bck=bck,
            cpu=cpu,
            quiet=quiet,
            fuzzy=fuzzy,
            fill=fill,
        )


    def _check_file_exists(self, path, description="file"):
        """
        Check if a file exists at the given path.

        Args:
            path (Path): The path to check.
            description (str): A description of the file for error messages.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        if not path.exists():
            raise FileNotFoundError(
                f"{description.capitalize()} '{str(path)}' does not exist."
            )


    def _validate_model(self, model, reference_img, remap):
        """
        Validate and set up the model, reference image, and remapping.

        Args:
            model (str or Path):
                Name of the predefined model or path to the custom model file.
            reference_img (str, optional):
                Path to a custom reference image.
            remap (str or dict, optional):
                Values used for remapping labels.
                It can be a string or a dictionary if it differs by side.

        Returns:
            Path to the model file.
            Path to the reference image.
            Code for remapping labels.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError:
                If a predefined model is not used
                and a reference image is not give.
        """
        PREDEFINED_MODELS = {
            "simple": {
                "weights": PROJECT_DIR / "lib" / "ensemble_hcvc.pth",
                "ref_img": PROJECT_DIR / "lib" / "ref_hcvc.mnc",
                "remap": { "left": {1:11, 2:12}, "right": {1:21, 2:22}}
            },
            "detailed": {
                "weights": PROJECT_DIR / "lib" / "ensemble_hcvc-ag.pth",
                "ref_img": PROJECT_DIR / "lib" / "ref_hcvc-ag.mnc",
                "remap": {
                    "left": {1:111, 2:112, 3:113, 4:121, 5:122, 6:123, 7:130},
                    "right": {1:211, 2:212, 3:213, 4:221, 5:222, 6:223, 7:230}
                }
            }
        }

        if model.lower() in PREDEFINED_MODELS:
            model_config = PREDEFINED_MODELS[model.lower()]
            model_path = model_config["weights"]
            ref_img_path = reference_img or model_config["ref_img"]
            remap = remap or model_config["remap"]

            self._check_file_exists(model_path, "Model file")
            self._check_file_exists(ref_img_path, "Reference image")
            return model_config["weights"], ref_img_path, remap

        model_path = Path(model)
        ref_img_path = Path(reference_img) if reference_img else None

        self._check_file_exists(model_path, "Model file")
        if ref_img_path:
            self._check_file_exists(ref_img_path, "Reference image")
        else:
            raise ValueError(
                "When using a non-predefined model, "
                "a path to a reference image is required."
            )

        return model_path, ref_img_path, remap


    #TODO: check about redundancy with self.params
    def _check_and_set_paths(self, in_img, output, work_dir):
        """Validate input/output/work_dir paths and set attributes."""
        updated_params = {}
        # Overwrite class attributes if arguments are not None
        if in_img:
            self.in_img = Path(in_img)
            updated_params["in_img"] = str(self.in_img)
        if output:
            self.output = Path(output)
            updated_params["output"] = str(self.output)
        if work_dir:
            self.work_dir = Path(work_dir)

        # Final checks
        if not self.in_img or not self.in_img.exists():
            raise FileNotFoundError(
                f"Input image '{str(self.in_img)}' does not exist."
            )

        if not self.output:
            self.output = in_img.parent
            warn(
                "Output was not set or found."
                f"'{str(self.output)}' will be used.",
                UserWarning
            )

        if self.output.is_dir():
            self._check_file_exists(self.output, "Output directory")
            self.output = self.output / f"{in_img.stem}_hcvc-segm.mnc"
            updated_params["output"] = str(self.output)

        if not self.work_dir or not self.work_dir.exists():
            self.work_dir = self.output.parent
            warn(
                "Working directory was not set or found. "
                f"'{str(self.work_dir)}' was set as working directory.",
                UserWarning
            )

        return updated_params


    def _process_side(self, minc, side, in_img, remap, clobber):
        """Helper to process one hemisphere."""
        side_img = self.work_dir / f"{in_img.stem}_{side}.mnc"
        side_seg = self.work_dir / f"{in_img.stem}_seg_{side}.mnc"
        side_rec = self.work_dir / f"{in_img.stem}_seg_rec_{side}.mnc"

        # Side transformation
        flip_xfm = self.work_dir / "flip.xfm"
        transform = str(flip_xfm) if side == "left" else None
        if not flip_xfm.exists() and transform:
            minc.param2xfm(str(flip_xfm), scales=[-1.,1.,1.])

        if side_img.exists() and clobber:
            side_img.unlink()

        minc.resample_smooth(
            str(in_img),
            str(side_img),
            like=str(self.ref_img),
            transform=transform,
        )

        if side_seg.exists() and clobber:
            side_seg.unlink()

        self.update_params(
            in_img=str(side_img),
            output=str(side_seg),
            verbose=False,
        )

        apply_cnn_segmentation(params=self.params)

        if side_rec.exists() and clobber:
            side_rec.unlink()

        remap_side = remap.get(side) if remap else None
        minc.resample_labels(
            str(side_seg),
            str(side_rec),
            like=str(in_img),
            transform=transform,
            remap=remap_side,
        )

        return side_rec

    def update_params(self, verbose=False, **kwargs):
        """
        Dynamically update the Params object with new parameters.

        Args:
            kwargs: Key-value pairs of parameters to update or add.
        """
        for key, value in kwargs.items():
            setattr(self.params, key, value)
        if verbose:
            print(f"\nUpdated params:\n{pformat(vars(self.params))}")
        self.has_run = False


    def run(
        self,
        in_img=None,
        output=None,
        work_dir=None,
        clobber=None,
        sep_by_sides=None
    ):
        """
        Run the CNN model using the stored parameters.

        Args:
            in_img (Path, optional): Path to the input image file.
            output (Path, optional): Path to save the output file.
        """
        clobber = clobber if clobber is not None else self.clobber
        sep_by_sides = (
            sep_by_sides
            if sep_by_sides is not None
            else self.sep_by_sides
        )

        updated_params = self._check_and_set_paths(in_img, output, work_dir)
        if updated_params:
            self.update_params(verbose=True, **updated_params)

        # Check if output exists, is a directory, or needs to be inferred
        if self.output.is_file():
            if clobber:
                warn(f"File {str(self.output)} already exists "
                     "and will be overwritten (clobber=True).\n",
                     UserWarning)
            else:
                self.has_run = True
                warn(f"File {str(self.output)} already exists.\n", UserWarning)
                return

        with mincTools() as minc:
            if sep_by_sides:
                # Left hemisphere
                left_rec = self._process_side(
                    minc=minc,
                    side="left",
                    in_img=self.in_img,
                    remap=self.remap,
                    clobber=clobber,
                )

                # Right hemisphere
                right_rec = self._process_side(
                    minc=minc,
                    side="right",
                    in_img=self.in_img,
                    remap=self.remap,
                    clobber=clobber,
                )

                # Combine Left & Right
                minc.calc(
                    inputs=[str(left_rec), str(right_rec)],
                    expression="A[0]>0?A[0]:A[1]",
                    output=str(self.output),
                    labels=True,
                )

                # Setback original in_img & output
                self.update_params(
                    verbose=False,
                    in_img=str(self.in_img),
                    output=str(self.output),
                )
            else:
                if self.remap:
                    seg = self.work_dir / f"{in_img.stem}_seg.mnc"
                    if seg.exists():
                        seg.unlink()
                    self.update_params(output=str(seg), verbose=False)
                    apply_cnn_segmentation(params=self.params)
                    self.update_params(output=str(self.output), verbose=False)
                    minc.resample_labels(
                        str(seg),
                        str(self.output),
                        like=str(self.in_img),
                        remap=self.remap
                    )
                else:
                    apply_cnn_segmentation(params=self.params)

        # Mark as done
        self.has_run = True


    def create_qc_img(self, output_qc=None, title=None, clobber=False):
        """
        Generate QC plot for a segmented image.

        Args:
            output_qc (str or Path, optional): Path to the output jpg file.
            title (str, optional): String to title the image.
            clobber (bool): Overwrite if existing.
        """
        # Validate that self.run has been completed
        if not self.has_run:
            raise RuntimeError("CNN must be run before generating QC.")

        title = title if title is not None else self.in_img.stem
        self.qc_img = (
            Path(output_qc) if output_qc
            else self.output.parent / f"{self.output.stem}.jpg"
        )


        if self.qc_img.exists():
            if clobber:
                self.qc_img.unlink()
            else:
                warn(
                    f"A QC image named '{self.qc_img.name}' already exists "
                    "and clobber option was not given.",
                    UserWarning
                )
                return

        lookup_table = PROJECT_DIR / "lib" / "labels.map"
        self._check_file_exists(lookup_table, "Lookup table for QC")

        remap_dict = {
            "simple": {11: 4, 12: 6, 21: 8, 22: 13},
            "detailed": {
                111: 1, 112: 2, 113: 3,
                121: 4, 122: 5, 123: 6,
                130: 7,
                211: 13, 212: 8, 213: 9,
                221: 14, 222: 11, 223: 12,
                230: 10,
            }
        }

        remap_codes = (
            remap_dict[self.model]
            if self.model in ("simple", "detailed") else None
        )

        print(f"remap_codes: {remap_codes}")

        with MincQC() as QC:
            QC.qc_function(
                input_img=self.in_img,
                input_labs=self.output,
                lookup_table=lookup_table,
                pik_name=self.qc_img,
                work_dir=self.work_dir,
                remap_dict=remap_codes,
                title=title,
            )
