#!/usr/bin/env python3
#
# Dr James Drewitt, 06/11/2020. last update: '19/07/2026'
# james.drewitt@bristol.ac.uk

import os
import sys
import subprocess
import argparse
from pathlib import Path

print('\n Go Magma');sys.stdout.flush()

project_dir = Path(__file__).resolve().parent
run_magma = project_dir / "MAGMA" / "run_magma.py"
launch_dir = Path.cwd()

DEFAULT_OPTIONS = {
    "contcar_file": "CONTCAR",
    "xyz_file": "movie.xyz",
    "xyz_numT": 1000,
    "T_step": 1,
    "alpha": "Si",
    "beta": "O",
    "r_cut_offs": [2.3],
    "save_detailed_analysis_data": False,
    "bond_angle_distribution": True,
    "partial_coordination": [4],
    "q_speciation": False,
    "carbonate_enabled": False,
    "carbonate_short_cutoff": 2.0,
    "carbonate_CC_cutoff": 2.0,
    "carbonate_long_min": 2.4,
    "carbonate_long_max": 2.9,
    "carbonate_angle_tolerance": 20.0,
    "carbonate_axial_tolerance": 20.0,
}


def _as_bool(value, line_number):
    if value.lower() in ("true", "yes", "1"):
        return True
    if value.lower() in ("false", "no", "0"):
        return False
    raise ValueError(f"line {line_number}: expected true or false")


def parse_text_config(config_text):
    """Parse a one-option-per-line MAGMA configuration."""
    options = DEFAULT_OPTIONS.copy()
    for line_number, line in enumerate(config_text.splitlines(), start=1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"line {line_number}: use option = value")
        name, value = (part.strip() for part in line.split("=", 1))
        if name not in DEFAULT_OPTIONS:
            raise ValueError(f"line {line_number}: unknown option '{name}'")
        if name in ("contcar_file", "xyz_file", "alpha", "beta"):
            options[name] = value
        elif name in ("xyz_numT", "T_step"):
            options[name] = int(value)
        elif name in ("r_cut_offs", "partial_coordination"):
            options[name] = [] if not value else [float(item) if name == "r_cut_offs" else int(item) for item in value.split(",")]
        elif name in ("save_detailed_analysis_data", "bond_angle_distribution", "q_speciation", "carbonate_enabled"):
            options[name] = _as_bool(value, line_number)
        else:
            options[name] = float(value)
    return options


def render_python_config(options):

    return f'''CONT_path = {("../" + options["contcar_file"])!r}
xyz = {("../" + options["xyz_file"])!r}
xyz_numT = {options["xyz_numT"]}
T_step = {options["T_step"]}
alpha = {options["alpha"]!r}
beta = {options["beta"]!r}
r_cut_offs = {options["r_cut_offs"]!r}
save_detailed_analysis_data = {int(options["save_detailed_analysis_data"])}
analysis = {{
    "bond_angle_distribution": {options["bond_angle_distribution"]!r},
    "partial_coordination": {options["partial_coordination"]!r},
    "q_speciation": {options["q_speciation"]!r},
    "carbonate": {{
        "enabled": {options["carbonate_enabled"]!r},
        "short_cutoff": {options["carbonate_short_cutoff"]},
        "cc_cutoff": {options["carbonate_CC_cutoff"]},
        "long_min": {options["carbonate_long_min"]},
        "long_max": {options["carbonate_long_max"]},
        "angle_tolerance": {options["carbonate_angle_tolerance"]},
        "axial_tolerance": {options["carbonate_axial_tolerance"]},
    }},
}}
'''

def default_text_config():
    def display(value):
        if isinstance(value, list):
            return ", ".join(map(str, value))
        return str(value).lower() if isinstance(value, bool) else str(value)
    return "\n".join(f"{name} = {display(value)}" for name, value in DEFAULT_OPTIONS.items()) + "\n"


def load_batch_input(input_path):
    """Read task folders and the optional shared configuration folder."""
    config_folder = input_path.parent
    tasks = []
    for line_number, raw_line in enumerate(input_path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.lower().startswith("config_folder"):
            name, separator, value = line.partition("=")
            if name.strip().lower() != "config_folder" or not separator or not value.strip():
                raise ValueError(f"{input_path}, line {line_number}: use config_folder = path")
            config_folder = Path(value.strip()).expanduser()
            if not config_folder.is_absolute():
                config_folder = input_path.parent / config_folder
        else:
            tasks.append(line)
    if not tasks:
        raise ValueError(f"{input_path} does not list any data folders")
    return tasks, config_folder

parser = argparse.ArgumentParser(description="Run MAGMA trajectory analysis.")

parser.add_argument(
    "--select",
    action="store_true",
    help="open a folder picker for serial runs",
)
parser.add_argument(
    "--confirm",
    action="store_true",
    help="show an editable configuration confirmation before the run (serial runs only)",
)
parser.add_argument(
    "-o",
    "--config-override",
    action="store_true",
    help="use the data folder's magma_config.txt if it exists, otherwise fall back to the shared batch config",
)
args = parser.parse_args()


def confirm_configuration(task_path, config_text, mode, number, total, has_saved_config):
    """Show an editable per-run configuration dialog"""
    import tkinter as tk
    from tkinter import messagebox

    result = {"config": None, "save": False}
    window = tk.Tk()
    window.title("MAGMA — confirm analysis")
    window.minsize(760, 600)

    tk.Label(
        window,
        text=f"{mode.title()} mode — run {number} of {total}\n{task_path}",
        justify="left",
        anchor="w",
    ).pack(fill="x", padx=12, pady=(12, 6))
    tk.Label(window, text="Edit this configuration if needed, then choose Run.", anchor="w").pack(
        fill="x", padx=12
    )

    frame = tk.Frame(window)
    frame.pack(fill="both", expand=True, padx=12, pady=8)
    scrollbar = tk.Scrollbar(frame)
    scrollbar.pack(side="right", fill="y")
    editor = tk.Text(frame, wrap="none", font=("Consolas", 10), yscrollcommand=scrollbar.set)
    editor.insert("1.0", config_text)
    editor.pack(side="left", fill="both", expand=True)
    scrollbar.config(command=editor.yview)

    save_changes = tk.BooleanVar(value=has_saved_config)
    tk.Checkbutton(
        window,
        text="Save these settings to this folder's magma_config.txt",
        variable=save_changes,
    ).pack(anchor="w", padx=12)

    buttons = tk.Frame(window)
    buttons.pack(fill="x", padx=12, pady=12)

    def run():
        candidate = editor.get("1.0", "end-1c")
        try:
            parse_text_config(candidate)
        except ValueError as error:
            messagebox.showerror("Invalid configuration", str(error))
            return
        result["config"] = candidate
        result["save"] = save_changes.get()
        window.destroy()

    tk.Button(buttons, text="Run", command=run, width=12).pack(side="right")
    tk.Button(buttons, text="Skip this folder", command=window.destroy, width=16).pack(side="right", padx=8)
    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.mainloop()
    return None if result["config"] is None else (result["config"], result["save"])

if args.select:
    from tkinter import Tk, filedialog

    picker = Tk()
    picker.withdraw()
    picker.attributes("-topmost", True)
    selected_folder = filedialog.askdirectory(
        title="Select the folder containing CONTCAR and out.xyz"
    )
    picker.destroy()
    if not selected_folder:
        sys.exit("No folder selected; MAGMA was not run.")
    tasks = [selected_folder]
else:
    try:
        tasks, batch_config_folder = load_batch_input(launch_dir / "input_magma.txt")
    except (OSError, ValueError) as error:
        sys.exit(str(error))

run_mode = "batch (sequential)" if len(tasks) > 1 else "serial"
print(f" Running in {run_mode} mode for {len(tasks)} folder(s).\n")
if args.confirm and run_mode != "serial":
    print(" --confirm is available only in serial mode; batch runs are non-interactive.\n")
if run_mode != "serial":
    config_mode = "independent configurations with shared fallback" if args.config_override else "shared configuration"
    print(f" Batch configuration mode: {config_mode}.\n")
    print(f" Shared configuration folder: {batch_config_folder}\n")

for task_number, task in enumerate(tasks, start=1):

        CWD = launch_dir
        task_path = Path(task).expanduser()
        if not task_path.is_absolute():
            task_path = Path(CWD) / task_path
        print(f"\n Running data folder {task_number} of {len(tasks)}: {task_path}\n")
        if run_mode != "serial" and not args.config_override:
            text_config_source = batch_config_folder / "magma_config.txt"
            if not text_config_source.is_file():
                sys.exit(
                    f"Configuration file not found: {text_config_source}\n"
                    "Set config_folder in input_magma.txt, create magma_config.txt there, "
                    "or use --config-override to use an independent configuration residing in the data-folder."
                )
        else:
            text_config_source = task_path / "magma_config.txt"
            if run_mode != "serial" and args.config_override and not text_config_source.is_file():
                shared_config_source = batch_config_folder / "magma_config.txt"
                if not shared_config_source.is_file():
                    sys.exit(
                        f"Independent configuration file not found: {text_config_source}\n"
                        f"Shared batch configuration file not found: {shared_config_source}"
                    )
                text_config_source = shared_config_source
                print(
                    f"Using shared batch configuration {text_config_source}.\n"
                )
        required_files = () if text_config_source.is_file() else ("CONTCAR", "movie.xyz")
        missing_inputs = [name for name in required_files if not (task_path / name).is_file()]
        if missing_inputs:
            sys.exit(f"{task_path} is missing required file(s): {', '.join(missing_inputs)}")
        taskdir = task_path / "MAGMA"

        # Create target directory if not exist
        if not os.path.exists(taskdir):
            os.makedirs(taskdir)
            print("\n Creating Directory " , taskdir, "\n");sys.stdout.flush()
        else:
            print("\n Creating Directory " , taskdir ,  " ... already exists\n");sys.stdout.flush()

        os.chdir(taskdir)#move into directory to run magma

        if text_config_source.is_file():
            public_config_text = text_config_source.read_text(encoding="utf-8")
            analysis_options = parse_text_config(public_config_text)
            config_text = render_python_config(analysis_options)
            print(f" Using independent configuration {text_config_source}\n")
        else:
            print(' Generating config file...\n' );sys.stdout.flush()
            public_config_text = default_text_config()
            analysis_options = parse_text_config(public_config_text)
            config_text = render_python_config(analysis_options)

        if args.confirm and run_mode == "serial":
            confirmed = confirm_configuration(
                task_path, public_config_text, run_mode, task_number, len(tasks), text_config_source.is_file()
            )
            if confirmed is None:
                print(f" Skipping {task_path}.\n")
                os.chdir(CWD)
                continue
            public_config_text, save_changes = confirmed
            analysis_options = parse_text_config(public_config_text)
            config_text = render_python_config(analysis_options)
            if save_changes:
                text_config_source.write_text(public_config_text, encoding="utf-8")

        child_environment = os.environ.copy()
        child_environment["MAGMA_CONFIG_TEXT"] = config_text
        result = subprocess.run(
            [sys.executable, str(run_magma)],
            cwd=taskdir,
            text=True,
            capture_output=True,
            env=child_environment,
        )
        if result.stdout:
            print(result.stdout, end="")
        if result.returncode:
            if result.stderr:
                print(result.stderr, file=sys.stderr, end="")
            sys.exit(result.returncode)

        os.chdir(CWD)#move back to original directory
