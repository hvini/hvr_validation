#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @author Sofia Fernandez-Lozano
# @date 20/11/2024
import argparse
import csv
from pathlib import Path
from tempfile import TemporaryDirectory
from tqdm import tqdm
from warnings import warn
from cnn.cnn_hcvc import CNN_HCVC
from minc.minc_calc_hvr import extract_hvr_from_img


def parse_csv(file_path, has_header, header_mapping=None):
    """
    Parse the input CSV file to extract subject information.

    Args:
        file_path (Path or str): Path to the CSV file.
        has_header (bool): Whether the CSV file has a header.
        header_mapping (dict): Custom header mapping for columns.

    Returns:
        List[Dict]: A list of dictionaries with columns as keys.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        # Try different heurestics of potential issues
        default_dir = Path("/app/data")
        #The parent dir of the filename was loaded to /app/data
        updated_path1 = default_dir / file_path.name

        #2 The complete path was loaded unto /app/data
        updated_path2 = default_dir / file_path

        if updated_path1.exists():
            file_path = updated_path1
        elif updated_path2.exists():
            file_path = updated_path2
        else:
            raise ValueError("CSV file was not found.")

    with file_path.open(mode='r') as csvfile:
        reader = csv.reader(csvfile)
        if has_header:
            headers = next(reader)
            # Dynamically map columns based on provided header_mapping
            img_path_col = headers.index(header_mapping["input_img_path"])
            try:
                subject_col = headers.index(header_mapping["subject_id"])
            except ValueError:
                subject_col = None
            try:
                session_col = headers.index(header_mapping["session_id"])
            except ValueError:
                session_col = None
            try:
                group_col = headers.index(header_mapping["group"])
            except ValueError:
                group_col = None
        else:
            row = next(reader)
            num_columns = len(row)

            # Mappings for variable column numbers
            column_mapping = {
                1: (None, None, 0, None),   # only image paths
                2: (0, None, 1, None),      # subj, paths
                3: (0, 1, 2, None),         # subj, sess, paths
                4: (0, 1, 3, 2)             # subj, sess, group, paths
            }

            # Assign columns based on the number of columns in the CSV
            if num_columns in column_mapping:
                subject_col, session_col, img_path_col, group_col = \
                    column_mapping[num_columns]
            else:
                raise ValueError(
                    f"CSV file has {num_columns} columns, "
                    "which is not supported."
                )

            # Return reader to the start of the file
            csvfile.seek(0)

        data = []
        for row in reader:
            subdict = {}
            if subject_col is not None:
                subdict["subject_id"] = row[subject_col]
            if session_col is not None:
                subdict["session_id"] = row[session_col]
            subdict["input_img_path"] = row[img_path_col]
            if group_col is not None:
                subdict["group"] = row[group_col]
            data.append(subdict)
        return data


def generate_output_paths(
    suffix="hcvc-segm",
    subject_id=None,
    session_id=None,
    input_img_path=None,
    segm_path=None,
    qc_path=None,
    qc_image=True,
):
    """
    Generate a special output name based on subject_id, session_id, and
    input_img_path.

    Args:
        subject_id (str, optional): Subject ID.
        session_id (str, optional): Session ID.
        input_img_path (Path, optional): Input image path.
        segm_path (Path, optional): Directory to save the output segmentations.
        qc_image (Bool): Whether to save QC images.
        qc_path (Path, optional): Directory to save the qc images.

    Returns:
        Dictionary: the different paths for a single subject
    """
    if not subject_id and not input_img_path:
        raise ValueError("Either subject_id or input_img_path must be given.")

    base_name = subject_id if subject_id else Path(input_img_path).stem

    if session_id:
        base_name = f"{base_name}_{session_id}"

    # Outputs
    outputs = {}
    segm_path = segm_path or Path("/app/output/segmentations")
    outputs["segmentation"] = segm_path / f"{base_name}_{suffix}.mnc"

    if qc_image:
        qc_path = qc_path or Path("/app/output/qc")
        outputs["qc_image"] = qc_path / f"{base_name}_{suffix}_qc.jpg"

    return outputs

def main(args):
    """
    Process a single or multiple subjects using the CNN_HCVC class.

    Args:
        args (Namespace): Parsed arguments.
    """
    # Work_dir set as argument otherwise create tmp directory
    if args.work_dir:
        work_dir = Path(args.work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        temp_dir = TemporaryDirectory()
        work_dir = Path(temp_dir.name)

    # Output_dir set as argument otherwise set default and give a warning.
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        # Use default directory inside the container.
        output_dir = Path("/app/output")
        warn("Output directory not specified. "
             "Defaulting to /app/output inside the container.\n"
             "Ensure you specify --output_dir if you want outputs to persist "
             "after the container stops.",
             UserWarning)

    try:
        if args.csv_file:
            header_mapping = {
                "subject_id": args.subject_id_col,
                "session_id": args.session_id_col,
                "input_img_path": args.input_img_col,
                "group": args.group_col
            }

            subjects = parse_csv(
                args.csv_file,
                args.has_header,
                header_mapping
            )
        else:
            # When no CSV, use the following heuristics:
            # Process the input_img_path argument.
            # If a CSV is not given, image path(s) is no longer optional.
            input_img_paths = [args.input_img_path] \
                if isinstance(args.input_img_path, str) \
                else args.input_img_path
            if input_img_paths is None:
                raise ValueError(
                    "A CSV file or a path to an input image has to be given."
                )

            # Number of paths is the number of processes that will be ran.
            num_images = len(input_img_paths)

            # Handle subject_id, session_id, and group arguments
            # Each can be either None, 1 (and be recycled),
            # or the same length as input_img_path
            subject_ids = args.subject_id.split(",") \
                if isinstance(args.subject_id, str) else args.subject_id
            if subject_ids is None:
                pass
            elif len(subject_ids) == 1:
                subject_ids = subject_ids * num_images
            elif len(subject_ids) != num_images:
                raise ValueError(
                    "Length of subject_id must be 1, None,"
                    " or match the length of input_img_path"
                    f" (got {len(subjects_ids)} subject_ids"
                    f" for {num_images} images)."
                )

            session_ids = args.session_id.split(",") \
                if isinstance(args.session_id, str) else args.session_id
            if session_ids is None:
                pass
            elif len(session_ids) == 1:
                session_ids = session_ids * num_images
            elif len(session_ids) != num_images:
                raise ValueError(
                    "Length of session_id must be 1, None,"
                    " or match the length of input_img_path"
                    f" (got {len(session_ids)} session_ids"
                    f" for {num_images} images)."
                )

            groups = args.group.split(",") \
                if isinstance(args.group, str) else args.group
            if groups is None:
                pass
            elif len(groups) == 1:
                groups = groups * num_images
            elif len(groups) != num_images:
                raise ValueError(
                    "Length of group must be 1, None,  or match the length of"
                    f" input_img_path (got {len(groups)} groups"
                    f" for {num_images} images)."
                )

            subjects = []
            for i in range(num_images):
                subject = {
                    "subject_id": subject_ids[i] \
                        if subject_ids is not None else None,
                    "session_id": session_ids[i] \
                        if session_ids is not None else None,
                    "input_img_path": input_img_paths[i],
                    "group": groups[i] if groups is not None else None
                }
                subjects.append(subject)

        # Initialize CNN_HCVC
        hcvc_segmenter = CNN_HCVC(
            model=args.model,
            work_dir=work_dir,
            clobber=args.clobber,
        )

        # Initialization of output directories and files
        seg_path = Path(args.seg_path) \
            if args.seg_path is not None \
            else output_dir / "segmentations"

        seg_path.mkdir(parents=True, exist_ok=True)

        qc_path = None
        if args.qc_images:
            qc_path = Path(args.qc_path) \
                if args.qc_path is not None \
                else output_dir / "qc"
            qc_path.mkdir(parents=True, exist_ok=True)

        if args.csv or args.hc_vols or args.vc_vols or args.hvr_csv:
            fieldnames = [
                "subject_id",
                "session_id",
                "group"
            ]

            for fieldname in fieldnames:
                if fieldname not in subjects[0]:
                    fieldnames.remove(fieldname)

            fieldnames.append("left")
            fieldnames.append("right")

            if args.csv or args.hvr_csv:
                # TODO: Make a function of this
                hvr_csv_path = Path(args.hvr_csv_path) \
                    if args.hvr_csv_path is not None \
                    else output_dir / "volumes" / "hvr_values.csv"

                hvr_csv_path.parent.mkdir(parents=True, exist_ok=True)

                if hvr_csv_path.exists():
                    warn(
                        f"The file {str(hvr_csv_path)} already exists. "
                        "It will be overwritten.",
                        UserWarning,
                    )
                    # TODO: Replace this with some other way: rename/append/idk
                    hvr_csv_path.unlink()

                with hvr_csv_path.open(mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(fieldnames)

            if args.csv or args.hc_vols:
                hc_csv_path = Path(args.hc_csv_path) \
                    if args.hc_csv_path is not None \
                    else output_dir / "volumes" / "hc_volumes.csv"

                hc_csv_path.parent.mkdir(parents=True, exist_ok=True)

                if hc_csv_path.exists():
                    warn(
                        f"The file {str(hc_csv_path)} already exists. "
                        "It will be overwritten.",
                        UserWarning,
                    )
                    # TODO: Replace this with some other way: rename/append/idk
                    hc_csv_path.unlink()

                with hc_csv_path.open(mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(fieldnames)

            if args.csv or args.vc_vols:
                vc_csv_path = Path(args.vc_csv_path) \
                    if args.vc_csv_path is not None \
                    else output_dir / "volumes" / "vc_volumes.csv"

                vc_csv_path.parent.mkdir(parents=True, exist_ok=True)

                if vc_csv_path.exists():
                    warn(
                        f"The file {str(vc_csv_path)} already exists. "
                        "It will be overwritten.",
                        UserWarning,
                    )
                    # TODO: Replace this with some other way: rename/append/idk
                    vc_csv_path.unlink()

                with vc_csv_path.open(mode='w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(fieldnames)

        # Process the subjects
        with tqdm(total=len(subjects), desc="Processing subjects") as pbar:
            for subject in subjects:
                input_img_path = Path(subject["input_img_path"])
                if not input_img_path.exists():
                    # Try different heurestics of potential issues
                    default_dir = Path("/app/data")
                    #The parent dir of the filename was loaded to /app/data
                    updated_path1 = default_dir / input_img_path.name

                    #2 The complete path was loaded unto /app/data
                    updated_path2 = default_dir / input_img_path

                    if updated_path1.exists():
                        input_img_path = updated_path1
                    elif updated_path2.exists():
                        input_img_path = updated_path2
                    else:
                        warn(
                            f"Input image {str(input_img_path)} does not exist. "
                            "Skipping.",
                            UserWarning
                        )
                        pbar.update(1)
                        continue

                subject_id = subject.get("subject_id")
                session_id = subject.get("session_id")
                group = subject.get("group")

                output_paths = generate_output_paths(
                    subject_id=subject_id,
                    session_id=session_id,
                    input_img_path=input_img_path,
                    segm_path=seg_path,
                    qc_image=args.qc_images,
                    qc_path=qc_path,
                )

                # Run the model
                hcvc_segmenter.run(
                    in_img=input_img_path,
                    output=output_paths["segmentation"]
                )

                # QC
                if args.qc_images:
                    # QC image title
                    title = subject_id
                    title = f"{title}_{session_id}" \
                        if title is not None and session_id is not None \
                        else None

                    hcvc_segmenter.create_qc_img(
                        output_qc=output_paths["qc_image"],
                        title=title,
                        clobber=args.clobber,
                    )

                # HC
                if args.csv or args.hc_vols or args.vc_vols or args.hvr_csv:
                    # TODO: Implement the option of it *not* being bilateral
                    # TODO: Implement the option of arguments for labels
                    vols_values = extract_hvr_from_img(
                        minc_path=output_paths["segmentation"],
                        # TODO: Double-check these two values with register
                        label_mapping={
                            "Left": {"HC": 11, "VC": 12},
                            "Right": {"HC": 21, "VC": 22},
                        },
                        bilateral=True,
                    )

                if args.csv or args.hc_vols:
                    row_csv = subject
                    row_csv["left"] = vols_values["Left"].get("HC")
                    row_csv["right"] =  vols_values["Right"].get("HC")

                    with hc_csv_path.open(mode='r+', newline='') as file:
                        reader = csv.reader(file)
                        header = next(reader, None)
                        row = [row_csv.get(column, "") for column in header]

                        writer = csv.writer(file)
                        file.seek(0, 2)
                        writer.writerow(row)

                if args.csv or args.vc_vols:
                    row_csv = subject
                    row_csv["left"] = vols_values["Left"].get("VC")
                    row_csv["right"] =  vols_values["Right"].get("VC")

                    with vc_csv_path.open(mode='r+', newline='') as file:
                        reader = csv.reader(file)
                        header = next(reader, None)
                        row = [row_csv.get(column, "") for column in header]

                        writer = csv.writer(file)
                        file.seek(0, 2)
                        writer.writerow(row)

                if args.csv or args.hvr_csv:
                    row_csv = subject
                    row_csv["left"] = vols_values["Left"].get("HVR")
                    row_csv["right"] =  vols_values["Right"].get("HVR")

                    with hvr_csv_path.open(mode='r+', newline='') as file:
                        reader = csv.reader(file)
                        header = next(reader, None)
                        row = [row_csv.get(column, "") for column in header]

                        writer = csv.writer(file)
                        file.seek(0, 2)
                        writer.writerow(row)

                pbar.update(1)
    finally:
        if not args.work_dir:
            temp_dir.cleanup()

def get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper script for CNN_HCVC processing.\n\n"
            "IMPORTANT:\n"
            "- All paths provided for input images, CSV files, and outputs "
            "must be located in bind-mounted directories.\n"
            "  This ensures data persists outside the container after it "
            "finishes running.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Input Options
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument(
        "-f", "--csv_file",
        help=(
            "Path to a CSV file with subject information. If no header is "
            "included, the columns will be guessed based on their count:\n"
            "- 1 column: Input image paths\n"
            "- 2 columns: Subject ID, Input image paths\n"
            "- 3 columns: Subject ID, Session ID, Input image paths\n"
            "- 4 columns: Subject ID, Session ID, Input image paths, Group\n"
            "IMPORTANT: This file must be located in a bind-mounted directory."
        )
    )
    input_group.add_argument(
        "-i", "--input_img_path",
        help=(
            "Path(s) to input image(s). "
            "Mandatory if no CSV file is given. It can be:\n"
            "- A single path\n"
            "- A comma-separated string of paths\n"
            "- A list of paths.\n"
            "IMPORTANT: All paths must lead to a bind-mounted directory."
        )
    )
    input_group.add_argument(
        "-s", "--subject_id",
        help=(
            "Optional. Subject ID(s). It can be:\n"
            "- A single string (recycled for all images)\n"
            "- A comma-separated string of IDs\n"
            "- A list of IDs.\n"
            "Must be 1 (to be recycled) or match the length of input_img_path."
        )
    )
    input_group.add_argument(
        "-v", "--session_id",
        help=(
            "Optional. Session ID(s). It can be:\n"
            "- A single string (recycled for all images)\n"
            "- A comma-separated string of IDs\n"
            "- A list of IDs.\n"
            "Must be 1 (to be recycled) or match the length of input_img_path."
        )
    )
    input_group.add_argument(
        "-g", "--group",
        help=(
            "Optional. Group(s). It can be:\n"
            "- A single string (recycled for all images)\n"
            "- A comma-separated string of groups\n"
            "- A list of groups.\n"
            "Must be 1 (to be recycled) or match the length of input_img_path."
        )
    )
    input_group.add_argument(
        "-m", "--model",
        default="simple",
        choices=["simple", "detailed"],
        help=(
            "Specify the CNN model:\n"
            "- 'simple': "
            "to segment complete Left/Right Hippocampus and Ventricles.\n"
            "- 'detailed': "
            "to segment Head/Body/Tail for HC & VC, as well as Amygdala."
        )
    )

    # Output Options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument(
        "-o", "--output_dir",
        help=(
            "Directory to save the output(s). "
            "This must be located in a bind-mounted directory.\n"
            "If missing all output will be saved on '/app/output' inside the "
            "container and might not persist after running."
        )
    )
    output_group.add_argument(
        "--vols", "--volumes_csv",
        action="store_true",
        dest="csv",
        help=("Save Hippocampal, Ventricular and HVR values into a CSV file.")
    )
    output_group.add_argument(
        "--hc_vols",
        action="store_true",
        help=("Save hippocampal volumes into a CSV file.")
    )
    output_group.add_argument(
        "--vc_vols",
        action="store_true",
        help=("Save ventricular volumes into a CSV file.")
    )
    output_group.add_argument(
        "--hvr", "--hvr_csv",
        action="store_true",
        dest="hvr_csv",
        help=("Calculate and save HVR values into a CSV file.")
    )
    output_group.add_argument(
        "--hc_path",
        dest="hc_csv_path",
        help=("Path to the CSV file with the HC values.")
    )
    output_group.add_argument(
        "--vc_path",
        dest="vc_csv_path",
        help=("Path to the CSV file with the VC values.")
    )
    output_group.add_argument(
        "--hvr_path",
        dest="hvr_csv_path",
        help=("Path to the CSV file with the HVR values.")
    )
    output_group.add_argument(
        "--qc", "--qc_images",
        action="store_true",
        dest="qc_images",
        help=("Include QC images of the segmentation(s) in the output.")
    )
    output_group.add_argument(
        "--qc_path",
        help=("Path to the directory where to store the QC images.")
    )
    output_group.add_argument(
        "--seg_path",
        help=("Path to the directory where to store the segmentations.")
    )
    output_group.add_argument(
        "-w", "--work_dir",
        help=(
            "Directory for intermediate files. "
            "If not provided, a temporary directory will be created.\n"
            "IMPORTANT: if specified, this must be a bind-mounted directory."
        )
    )
    output_group.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite existing files if they already exist."
    )

    # CSV Column Options
    csv_group = parser.add_argument_group("CSV Column Options")
    csv_group.add_argument(
        "--has_header",
        action="store_true",
        help="Indicates if the CSV file has a header row."
    )
    csv_group.add_argument(
        "--subject_id_col",
        default="subject_id",
        help=(
            "Column name for subject IDs in the CSV file. "
            "Defaults to 'subject_id'."
        )
    )
    csv_group.add_argument(
        "--session_id_col",
        default="session_id",
        help=(
            "Column name for session IDs in the CSV file. "
            "Defaults to 'session_id'."
        )
    )
    csv_group.add_argument(
        "--input_img_col",
        default="input_img_path",
        help=(
            "Column name for input image paths in the CSV file. "
            "Defaults to 'input_img_path'."
        )
    )
    csv_group.add_argument(
        "--group_col",
        default="group",
        help=(
            "Column name for group variable in the CSV file. "
            "Defaults to 'group'."
        )
    )

    return parser


def check_required_args(args):
    if not args.csv_file and not args.input_img_path:
        raise ValueError(
            "You must specify either --csv_file or --input_img_path."
        )
    if args.csv_file and args.input_img_path:
        raise ValueError(
            "Specify only one of --csv_file or --input_img_path, not both."
        )
    if args.csv_file and (args.subject_id or args.session_id or args.group):
        warn(
            "You specified both --csv_file and --subject_id, --session_id, "
            "and/or --group.\n"
            "When using a CSV file, the latter options are ignored. "
            "Set metadata using the CSV Column Options. See --help for info.",
            UserWarning
        )

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    check_required_args(args)
    main(args)
