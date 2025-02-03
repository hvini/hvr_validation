# -*- coding: utf-8 -*-

# @author Sofia Fernandez
# @date 2024/11/22

# Adaptation from older HC QC Perl script

import argparse
from tempfile import TemporaryDirectory
from pathlib import Path
from random import randint
from minc.minc_tools import mincTools

class MincQC(mincTools):
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


    def resample_files(
        self,
        in_mri,
        out_mri,
        in_labels=None,
        out_labels=None,
        remap_dict=None
    ):
        """
        Resample MRI and label files to focus on Mediotemporal area.

        Args:
            in_mri (Path): The path to the input MRI.
            out_mri (Path): The path where to save the resampled MRI.
            in_labels (Path; optional): The path to the input labels.
            out_labels (Path; optional): The path where to save the
                resampled labels.
            remap_dict (dict; optional): Dictionary with remapping values for
                the labels
        """
        # Resample input image
        resample_mri_cmd = [
            "mincresample",
            str(in_mri), str(out_mri),
            "-trilinear",
            "-q",
            "-fill",
            "-fillvalue", str(0.),
            "-step", str(1.), str(1.), str(1.),
            "-start", str(-65.), str(-64.), str(-55.),
            "-nelements", str(131), str(98), str(84),
        ]

        self.command(
            resample_mri_cmd,
            inputs=[str(in_mri)],
            outputs=[str(out_mri)],
            verbose=self.verbose
        )

        # Resample labels
        if in_labels and out_labels:
            # TODO: Move this to outside this script
            self.resample_labels(
                input=str(in_labels),
                output=str(out_labels),
                like=str(out_mri),
                datatype="byte",
                remap=remap_dict,
            )


    def normalize_mri(self, in_mri, out_mri):
        """
        Normalize MRI intensities to 0-100.

        Args:
            in_mri (Path): The path to the input MRI.
            out_mri (Path): The path where to save the normalized MRI.
        """

        # Command
        norm_mri_cmd = [
            "mincnorm",
            "-cutoff", str(0.5),
            str(in_mri),
            str(out_mri)
        ]

        self.command(
            norm_mri_cmd,
            inputs=[str(in_mri)],
            outputs=[str(out_mri)],
            verbose=self.verbose
        )


    def final_lookup(
            self,
            in_mri,
            in_labels,
            mapping,
            work_dir,
            out_img
    ):
        """
        Apply final mapping and create a single image for QC.

        Args:
            in_mri (Path): The path to the input MRI.
            in_labels (Path): The path to the input labels.
            mapping (Path): The path to the mapping file used for lookup.
            work_dir (Path): The path to the working directory where to save
                intermediate files.
            out_img (Path): The path where to save the QC MRI.
        """
        prefix = out_img.stem

        # Grey lookup MRI
        grey_mri = Path(work_dir) / f"{prefix}_grey_mri.mnc"
        grey_cmd = [
            "minclookup",
            "-clobber",
            "-grey",
            "-range", str(20), str(90),
            str(in_mri),
            str(grey_mri),
        ]

        # Labels lookup with mapping
        final_labs = Path(work_dir) / f"{prefix}_final_labs.mnc"
        labs_cmd = [
            "minclookup",
            "-clobber",
            "-discrete",
            "-lookup_table", str(mapping),
            str(in_labels),
            str(final_labs),
        ]

        # Max of both images
        merge_cmd = [
            "mincmath",
            "-nocheck_dimensions",
            "-max",
            str(grey_mri),
            str(final_labs),
            str(out_img),
        ]

        # Running commands
        self.command(
            grey_cmd,
            inputs=[str(in_mri)],
            outputs=[str(grey_mri)],
            verbose=self.verbose
        )

        self.command(
            labs_cmd,
            inputs=[str(in_labels), str(mapping)],
            outputs=[str(final_labs)],
            verbose=self.verbose
        )

        self.command(
            merge_cmd,
            inputs=[str(grey_mri), str(final_labs)],
            outputs=[str(out_img)],
            verbose=self.verbose
        )


    def minc_pik(self, qc_img, work_dir, output, text=None):
        """
        Create a layout of slices for QC. Specific to medio-temporal area.

        Args:
            qc_img (Path): The path to the input MINC file.
            work_dir (Path): The path to the directory
                used for intermediate files.
            output (Path): The path to the final output image file.
            text (str; optional): Title text to add to the output image.
        """
        prefix = output.stem

        slices = {
            "transverse": [26, 30, 34, 38, 42, 44, 48],
            "sagittal1": [28, 31, 34, 37, 40, 43, 47],
            "sagittal2": [102, 99, 96, 93, 90, 87, 84],
            "coronal": [55, 50, 45, 40, 35, 30, 25],
        }

        miff_files = []  # All intermediate .miff files for montage
        small_tile_size = 200  # Tile size for the montage

        # Loop through slice types and generate commands
        for _orientation, slice_list in slices.items():
            for slice_number in slice_list:
                # Set orientation-specific prefix
                _prefix = _orientation[0].upper()
                orientation = (
                    _orientation if _prefix != "S" else _orientation[:-1]
                )
                _prefix = f"{prefix}_{_prefix}"
                miff_file = work_dir / f"{_prefix}{slice_number}.miff"

                # Build the command for mincpik
                mincpik_cmd = [
                    "mincpik",
                    "-scale", "1",
                    f"-{orientation}",
                    "-slice",
                    str(slice_number),
                    "-clobber",
                    str(qc_img),
                    str(miff_file)
                ]

                # Run the command
                self.command(
                    mincpik_cmd,
                    inputs=[str(qc_img)],
                    outputs=[str(miff_file)],
                    verbose=self.verbose
                )

                # Add the resulting miff file to the montage command list
                miff_files.append(str(miff_file))

        # Create montage from the generated miff files
        _output = work_dir / f"{prefix}_mont.miff"
        montage_cmd = [
            "montage",
            "-tile",
            "7x4",
            "-background",
            "grey10",
            "-geometry",
            f"{small_tile_size}x{small_tile_size}+1+1"
        ] + miff_files + [
            str(_output)
        ]

        # Run the command
        self.command(
            montage_cmd,
            inputs=miff_files,
            outputs=[str(_output)],
            verbose=self.verbose
        )

        # Add optional title text to the montage
        if text:
            draw_cmd = [
                "convert",
                str(_output),
                "-box",
                "white",
                "-draw",
                f'text 2,15 "{text}"',
                str(output)
            ]
        else:
            draw_cmd = [
                "convert",
                str(_output),
                str(output)
            ]

        # Run the command
        self.command(
            draw_cmd,
            inputs=[str(_output)],
            outputs=[str(output)],
            verbose=self.verbose
        )


    def qc_function(
        self,
        input_img,
        input_labs,
        lookup_table,
        pik_name,
        title=None,
        work_dir=None,
        remap_dict=None,
    ):
        """
        Performs QC by reshaping, looking up intensities, and picking slices.

        Args:
            input_img (str, Path): Path to input MINC brain image
            input_labs (str, Path): Path to input MINC segmentation labels
            pik_name (str, Path): Path for the pik output file
            title (str; optional): Title for the qc image.
            lookup_table (str, Path): Path for the lookup file
            work_dir (str, Path; optional):
                Path to directory to save intermediate files
            remap_dict (dict; optional):
                Remapping values for re-labelling before lookup
        """
        # Consolidate and check paths
        input_img = Path(input_img)
        self._check_file_exists(input_img, description="input image")

        input_labs = Path(input_labs)
        self._check_file_exists(input_labs, description="input labels")

        if work_dir is not None:
            work_dir = Path(work_dir)
            self._check_file_exists(work_dir, description="working directory")
            temp_dir = None
        else:
            temp_dir = TemporaryDirectory()
            work_dir = Path(temp_dir.name)


        lookup_table = Path(lookup_table)
        self._check_file_exists(lookup_table, description="lookup table")

        pik_name = Path(pik_name)
        self._check_file_exists(
            pik_name.parent,
            description="output directory"
        )

        _title = (
            f"{title}_"
            if title is not None
            else "".join(str(randint(0, 9)) for _ in range(10))
        )

        try:
            # Resample files
            resampled_img = work_dir / f"{_title}resampled_img.mnc"
            resampled_labs = work_dir / f"{_title}resampled_lab.mnc"
            self.resample_files(
                in_mri=input_img,
                in_labels=input_labs,
                out_mri=resampled_img,
                out_labels=resampled_labs,
                remap_dict=remap_dict,
            )

            # Normalize MRI
            normalized_img = work_dir / f"{_title}normalized_img.mnc"
            self.normalize_mri(
                in_mri=resampled_img,
                out_mri=normalized_img,
            )

            # QC MINC
            qc_img = work_dir / f"{_title}qc.mnc"
            self.final_lookup(
                in_mri=normalized_img,
                in_labels=resampled_labs,
                mapping=lookup_table,
                work_dir=work_dir,
                out_img=qc_img,
            )

            # MINCPIK
            self.minc_pik(
                qc_img=qc_img,
                work_dir=work_dir,
                output=pik_name,
                text=title,
            )

        finally:
            if temp_dir is not None:
                temp_dir.cleanup()


def main():
    parser = argparse.ArgumentParser(description="QC function for MINC files.")
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument(
        "input_img",
        help="Path to the input MINC file."
    )
    input_group.add_argument(
        "-l", "--labels",
        dest="input_labs",
        help="Path to the input MINC segmentation labels."
    )
    input_group.add_argument(
        "--lut", "--lookup-table",
        dest="lookup_table",
        help="Path for the lookup table file."
    )
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument(
        "-o", "--output",
        dest="pik_name",
        help="Path and filename for the output QC image."
    )
    optional_group = parser.add_argument_group("Optional Options")
    optional_group.add_argument(
        "-w", "--work_dir",
        help=(
            "Directory for intermediate files. "
            "If not provided, a temporary directory will be created."
        )
    )
    optional_group.add_argument(
        "-m", "--remap_dict",
        help="Remapping values for re-labelling before lookup."
    )

    # Create an instance of MincQC and run the QC function
    qc_tool = MincQC()
    qc_tool.qc_function(
        input_img=args.input_img,
        input_labs=args.input_labs,
        lookup_table=args.lookup_table,
        pik_name=args.pik_name,
        work_dir=args.work_dir,
        remap_dict=args.remap_dict,
    )

if __name__ == "__main__":
    main()
