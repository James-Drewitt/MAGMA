#!/usr/bin/env python3
#
# Dr James Drewitt, 19/07/2026. Last update: "19/07/2026"
#
#Calculate population and lifetime of carbonate structural units.

from itertools import combinations
from pathlib import Path
import os

import numpy as np


def _minimum_image(vectors, box_length):
    return vectors - box_length * np.rint(vectors / box_length)


def _lifetimes(frame_units, step):
    #Return consecutive lifetimes for labelled units in frames.
    active = {}
    runs = []
    for units in frame_units:
        for unit in tuple(active):
            if unit not in units:
                runs.append(active.pop(unit))
        for unit in units:
            active[unit] = active.get(unit, 0) + 1
    # Units still active at the trajectory boundary are right-censored, but
    # their observed duration belongs in the reported distribution.
    runs.extend(active.values())
    return [run * step for run in runs]


def _write_lifetime(path, unit_name, lifetimes):
    if lifetimes:
        values = [
            f"mean lifetime of {unit_name} = {np.mean(lifetimes)} timesteps",
            f"median lifetime = {np.median(lifetimes)} timesteps",
            f"min lifetime = {min(lifetimes)} timesteps",
            f"max lifetime = {max(lifetimes)} timesteps",
        ]
    else:
        values = [f"no {unit_name} units were found"]
    np.savetxt(path, np.asarray(values, dtype=object), fmt="%s")


def _append_cn_summary(output_dir, prefix, CO3plus1_counts, CO4_counts, C2O5_counts, carbon_counts):
    """Append carbonate populations to the coordination-number average file."""
    total_CO4 = sum(CO4_counts)
    total_CO3plus1 = sum(CO3plus1_counts)
    total_carbons = sum(carbon_counts)
    total_C2O5 = sum(C2O5_counts)
    CO3plus1_percent = 100.0 * total_CO3plus1 / total_CO4 if total_CO4 else 0.0
    mean_C2O5 = float(np.mean(C2O5_counts)) if C2O5_counts else 0.0
    C2O5_carbon_percent = 200.0 * total_C2O5 / total_carbons if total_carbons else 0.0

    with (output_dir / f"{prefix}-CN-av.dat").open("a", encoding="utf-8") as handle:
        handle.write("\ncarbonate_summary\n")
        handle.write("CO3+1 (trigonal planar CO3 unit with an additional long axial C--O bond)\n")
        handle.write(f"CO3+1 fraction of all CO4: {CO3plus1_percent:.4f} %\n\n")
        handle.write(f"mean number of C2O5 units: {mean_C2O5:.6f}\n")
        handle.write(f"carbon atoms in C2O5 units: {C2O5_carbon_percent:.4f} %\n")


def carbonate_units(
    filename, t_step, alpha, beta, box_length, xyz_num_t, working_dir,
    short_cutoff=2.0, CC_cutoff=2.0, long_min=2.4, long_max=2.9,
    angle_tolerance=20.0, axial_tolerance=20.0,
):
    #Find CO3+1 and C2O5 units in an XYZ trajectory.

    #CO3+1 is a trigonal-planar carbonate unit with three short C--O bonds and 
    #one long approximately axial C--O bond.

    #C2O5 is defined as either a pair of planar CO3 groups sharing one short 
    #C--O--C bridge, or a short homopolar C--C pair with three and two distinct 
    #short C--O neighbours (O3--C=C--O2).
    
    if (alpha, beta) != ("C", "O"):
        raise ValueError("Carbonate-unit analysis requires alpha='C' and beta='O'.")
    if not 0 < short_cutoff < long_min <= long_max:
        raise ValueError("Carbonate cut-offs must satisfy 0 < short < long_min <= long_max.")
    if CC_cutoff <= 0:
        raise ValueError("The carbonate C--C cutoff must be positive.")

    with open(filename, encoding="utf-8") as handle:
        lines = handle.readlines()
    n_atoms = int(lines[0])
    frame_size = n_atoms + 2
    n_traj = len(lines) // frame_size
    if n_traj > xyz_num_t:
        lines = lines[(n_traj - xyz_num_t) * frame_size:]
        n_traj = xyz_num_t

    sin_angle = np.sin(np.deg2rad(angle_tolerance))
    cos_axial = np.cos(np.deg2rad(axial_tolerance))
    CO3plus1_frames = []
    C2O5_frames = []
    C2O5_bridge_frames = []
    C2O5_homopolar_frames = []
    CO3plus1_rows = []
    C2O5_rows = []
    CO3plus1_counts = []
    CO4_counts = []
    C2O5_counts = []
    carbon_counts = []

    for frame in range(0, n_traj, t_step):
        atom_rows = [line.split() for line in lines[2 + frame * frame_size:(frame + 1) * frame_size]]
        carbons = [(f"C{i + 1}", *map(float, row[1:])) for i, row in enumerate(atom_rows) if row[0] == alpha]
        oxygens = [(f"O{i + 1}", *map(float, row[1:])) for i, row in enumerate(atom_rows) if row[0] == beta]
        if not carbons or not oxygens:
            raise ValueError(f"Frame {frame + 1} does not contain both C and O atoms.")

        C_xyz = np.asarray([row[1:] for row in carbons], dtype=float)
        O_xyz = np.asarray([row[1:] for row in oxygens], dtype=float)
        vectors = _minimum_image(O_xyz[None, :, :] - C_xyz[:, None, :], box_length)
        distances = np.linalg.norm(vectors, axis=2)
        short_oxygen_labels = {
            c_index: frozenset(oxygens[index][0] for index in np.flatnonzero(distances[c_index] <= short_cutoff))
            for c_index in range(len(carbons))
        }

        trigonal = {}
        CO3plus1 = set()
        CO4_count = 0
        for c_index, carbon in enumerate(carbons):
            short = np.flatnonzero(distances[c_index] <= short_cutoff)
            long_bonds = np.flatnonzero((distances[c_index] >= long_min) & (distances[c_index] <= long_max))
            if len(short) == 3 and len(long_bonds) == 1:
                CO4_count += 1
            if len(short) != 3:
                continue
            short_vectors = vectors[c_index, short]
            norms = np.linalg.norm(short_vectors, axis=1)
            unit_vectors = short_vectors / norms[:, None]
            pair_angles = [
                np.rad2deg(np.arccos(np.clip(np.dot(unit_vectors[a], unit_vectors[b]), -1.0, 1.0)))
                for a, b in combinations(range(3), 2)
            ]
            normal = np.cross(unit_vectors[0], unit_vectors[1])
            normal_norm = np.linalg.norm(normal)
            if normal_norm == 0 or any(abs(angle - 120.0) > angle_tolerance for angle in pair_angles):
                continue
            normal /= normal_norm
            if abs(np.dot(normal, unit_vectors[2])) > sin_angle:
                continue

            short_labels = frozenset(oxygens[index][0] for index in short)
            trigonal[c_index] = short_labels
            axial = [
                index for index in long_bonds
                if abs(np.dot(vectors[c_index, index] / distances[c_index, index], normal)) >= cos_axial
            ]
            if len(long_bonds) == 1 and len(axial) == 1:
                CO3plus1.add(f"{carbon[0]}-{oxygens[axial[0]][0]}")

        C2O5_bridge = set()
        for first, second in combinations(trigonal, 2):
            shared = trigonal[first] & trigonal[second]
            if len(shared) == 1:
                bridge = next(iter(shared))
                C2O5_bridge.add(f"{carbons[first][0]}-{bridge}-{carbons[second][0]}")

        C2O5_homopolar = set()
        CC_vectors = _minimum_image(C_xyz[None, :, :] - C_xyz[:, None, :], box_length)
        CC_distances = np.linalg.norm(CC_vectors, axis=2)
        for first, second in combinations(range(len(carbons)), 2):
            first_oxygens = short_oxygen_labels[first]
            second_oxygens = short_oxygen_labels[second]
            if (
                CC_distances[first, second] <= CC_cutoff
                and sorted((len(first_oxygens), len(second_oxygens))) == [2, 3]
                and len(first_oxygens | second_oxygens) == 5
            ):
                C2O5_homopolar.add(f"{carbons[first][0]}={carbons[second][0]}")

        C2O5 = (
            {f"oxygen_bridge:{unit}" for unit in C2O5_bridge}
            | {f"homopolar_C=C:{unit}" for unit in C2O5_homopolar}
        )

        frame_number = frame + 1
        CO3plus1_frames.append(CO3plus1)
        C2O5_frames.append(C2O5)
        C2O5_bridge_frames.append(C2O5_bridge)
        C2O5_homopolar_frames.append(C2O5_homopolar)
        CO3plus1_counts.append([frame_number, len(CO3plus1)])
        CO4_counts.append(CO4_count)
        C2O5_counts.append([frame_number, len(C2O5), len(C2O5_bridge), len(C2O5_homopolar)])
        carbon_counts.append(len(carbons))
        CO3plus1_rows.extend([[frame_number, unit] for unit in sorted(CO3plus1)])
        C2O5_rows.extend([[frame_number, "oxygen_bridge", unit] for unit in sorted(C2O5_bridge)])
        C2O5_rows.extend([[frame_number, "homopolar_C=C", unit] for unit in sorted(C2O5_homopolar)])

    output_dir = Path(os.getcwd()) / working_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{alpha}-{beta}"
    np.savetxt(
        output_dir / f"{prefix}_CO3plus1_counts.dat",
        np.asarray([["trajectory", "CO3plus1_count"], *CO3plus1_counts], dtype=object), fmt="%s",
    )
    np.savetxt(
        output_dir / f"{prefix}_CO3plus1_membership.dat",
        np.asarray([["trajectory", "unit_identifier"], *CO3plus1_rows], dtype=object), fmt="%s",
    )
    np.savetxt(
        output_dir / f"{prefix}_C2O5_counts.dat",
        np.asarray(
            [["trajectory", "C2O5_count", "oxygen_bridge_count", "homopolar_C=C_count"], *C2O5_counts],
            dtype=object,
        ), fmt="%s",
    )
    np.savetxt(
        output_dir / f"{prefix}_C2O5_membership.dat",
        np.asarray([["trajectory", "origin", "unit_identifier"], *C2O5_rows], dtype=object), fmt="%s",
    )
    _write_lifetime(output_dir / f"{prefix}_CO3plus1_lifetime.dat", "CO3+1", _lifetimes(CO3plus1_frames, t_step))
    _write_lifetime(output_dir / f"{prefix}_C2O5_lifetime.dat", "C2O5", _lifetimes(C2O5_frames, t_step))
    _write_lifetime(
        output_dir / f"{prefix}_C2O5_oxygen_bridge_lifetime.dat",
        "oxygen-bridge C2O5", _lifetimes(C2O5_bridge_frames, t_step),
    )
    _write_lifetime(
        output_dir / f"{prefix}_C2O5_homopolar_CC_lifetime.dat",
        "homopolar C2O5", _lifetimes(C2O5_homopolar_frames, t_step),
    )
    _append_cn_summary(
        output_dir, prefix,
        [count for _, count in CO3plus1_counts], CO4_counts,
        [count for _, count, _, _ in C2O5_counts], carbon_counts,
    )
    return CO3plus1_frames, C2O5_frames
