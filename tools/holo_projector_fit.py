#!/usr/bin/env python3
"""Holo projector poses and meshes for the dome model.

Two jobs, both offline (numpy only):

  poses   Find the three holo projector holes in meshes/dome/base_dome.obj and write
          meshes/dome/holo_constants.xacro: for each projector the pivot point and the
          orientation of its frame, in the dome link's frame. The cowl's flange is inside
          the dome, flat against the inner skin, and its cone comes out through the hole.
          So the lens axis is the skin's normal at the hole (the holes are not cut quite
          square to the skin), and the pivot (the ball's centre) is under the middle of
          the hole, level with the flange's base. Re-run after base_dome.obj is
          re-exported.

  meshes  Re-origin the two CSL HoloEye STLs: the body so its origin is the centre of its
          ball and the cowl so its origin is the centre of its flange's base, each with the
          lens axis along +Z. Units stay millimetres (scale 0.001 in the URDF).

The projector frame: +Z out of the dome along the skin's normal, +X horizontal (to the
viewer's right when looking into the lens), +Y towards the top of the dome.

    python3 tools/holo_projector_fit.py poses
    python3 tools/holo_projector_fit.py meshes BODY.stl COWL.stl
"""

import argparse
import math
import pathlib
import struct

import numpy as np

PACKAGE = pathlib.Path(__file__).resolve().parent.parent
DOME_DIR = PACKAGE / 'meshes' / 'dome'
DOME_SCALE = 0.01          # base_dome.obj is in centimetres
DOME_UP = np.array([0.0, 1.0, 0.0])

# Roughly where each hole is, in base_dome.obj units. Only used to pick the hole; the
# pose comes from the fit.
SEEDS = {
    'hp_front': (-20.2, 6.7, 9.6),
    'hp_top': (8.2, 24.4, 5.3),
    'hp_rear': (21.9, 7.4, -3.7),
}
HOLE_RADIUS = (2.9, 3.6)   # accepted distance of a hole-wall face from the hole's axis
HOLE_DEPTH = 1.5           # and along it
HOLE_FIT_TOLERANCE = 0.1   # worst distance of a wall vertex from the fitted cylinder


def fit_sphere(points):
    m = np.c_[2 * points, np.ones(len(points))]
    x = np.linalg.lstsq(m, (points ** 2).sum(1), rcond=None)[0]
    centre = x[:3]
    return centre, math.sqrt(x[3] + centre @ centre)


def load_obj(path):
    vertices, faces = [], []
    for line in open(path):
        if line.startswith('v '):
            vertices.append([float(v) for v in line.split()[1:4]])
        elif line.startswith('f '):
            idx = [int(t.split('/')[0]) - 1 for t in line.split()[1:]]
            faces += [[idx[0], idx[i], idx[i + 1]] for i in range(1, len(idx) - 1)]
    return np.array(vertices)[np.array(faces)]


def fit_plane(points):
    """Return (centroid, unit normal) of the plane through points."""
    centroid = points.mean(0)
    return centroid, np.linalg.svd(points - centroid)[2][2]


def fit_hole(name, tris, centres, normals, areas, seed, dome_centre, flange):
    """Return (pivot, normal, radius, gap) for the round hole near seed.

    normal is the dome skin's outward normal at the hole and pivot the ball's centre;
    flange is the cowl flange's (radius, thickness) in the dome mesh's units; gap is how
    far the flange's top is from the inner skin at the hole's edge.
    """
    out = (seed - dome_centre) / np.linalg.norm(seed - dome_centre)
    centre = seed
    # The fit picks the wall faces by distance from the hole's axis, so repeat it from its
    # own answer until a rough seed no longer matters
    for _ in range(4):
        d = centres - centre
        h = d @ out
        rad = np.linalg.norm(d - h[:, None] * out, axis=1)
        wall = ((areas > 1e-9) & (rad > HOLE_RADIUS[0]) & (rad < HOLE_RADIUS[1])
                & (np.abs(h) < HOLE_DEPTH) & (np.abs(normals @ out) < 0.3))
        if wall.sum() < 50:
            raise SystemExit(f'{name}: no hole found near {np.round(centre, 2)}')

        # The cut's axis is the direction every wall normal is perpendicular to
        cut = np.linalg.svd(normals[wall] * np.sqrt(areas[wall])[:, None])[2][2]
        cut *= np.sign(cut @ out)

        # Circle fit of the wall vertices in the plane across the cut
        points = np.unique(tris[wall].reshape(-1, 3), axis=0)
        e1 = np.cross(cut, DOME_UP)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(cut, e1)
        xy = np.c_[(points - centre) @ e1, (points - centre) @ e2]
        x = np.linalg.lstsq(np.c_[2 * xy, np.ones(len(xy))], (xy ** 2).sum(1),
                            rcond=None)[0]
        radius = math.sqrt(x[2] + x[0] ** 2 + x[1] ** 2)
        residual = np.abs(np.linalg.norm(xy - x[:2], axis=1) - radius).max()
        coverage = len(set(np.floor(np.arctan2(xy[:, 1] - x[1], xy[:, 0] - x[0])
                                    / (math.pi / 6)).astype(int)))
        centre = centre + x[0] * e1 + x[1] * e2
        out = cut

    if (residual > HOLE_FIT_TOLERANCE or coverage < 12
            or not HOLE_RADIUS[0] < radius < HOLE_RADIUS[1]):
        raise SystemExit(f'{name}: the faces near {np.round(seed, 2)} are not a round hole '
                         f'(radius {radius:.3f}, worst vertex {residual:.3f} off, '
                         f'{coverage}/12 of the way round); fix its entry in SEEDS')

    # The hole is not cut square to the skin, and the flange lies against the skin. So the
    # projector's axis is the skin's normal: the normal of the plane through the hole's
    # inner rim. The wall's vertices are the two rims; the inner one is nearer the centre.
    from_centre = np.linalg.norm(points - dome_centre, axis=1)
    inner_rim = points[from_centre < (from_centre.min() + from_centre.max()) / 2]
    rim_centre, normal = fit_plane(inner_rim)
    normal *= np.sign(normal @ cut)
    hole_centre = centre + ((rim_centre - centre) @ normal) / (cut @ normal) * cut

    # A flat flange inside a sphere touches at its edge, so its top is below the inner
    # rim by the skin's fall between the hole's radius and the flange's
    flange_radius, flange_thickness = flange
    skin_radius = np.linalg.norm(inner_rim - dome_centre, axis=1).mean()
    gap = (math.sqrt(skin_radius ** 2 - radius ** 2)
           - math.sqrt(skin_radius ** 2 - flange_radius ** 2))

    # Pivot: under the middle of the hole, level with the flange's base
    return hole_centre - (gap + flange_thickness) * normal, normal, radius, gap


def cowl_flange():
    """(radius, thickness) of the cowl's flange, in the dome mesh's units."""
    points = read_stl(DOME_DIR / 'holo_projector_cowl.stl')['v'].reshape(-1, 3)
    r = np.linalg.norm(points[:, :2], axis=1)
    mm = 0.001 / DOME_SCALE
    return r.max() * mm, points[r > 0.9 * r.max(), 2].max() * mm


def frame_rpy(axis):
    """URDF rpy of the frame with +Z along axis, +X horizontal, +Y up the dome."""
    x = np.cross(DOME_UP, axis)
    if np.linalg.norm(x) < 1e-3:
        raise SystemExit('a hole pointing straight up has no horizontal +X; choose one here')
    x /= np.linalg.norm(x)
    y = np.cross(axis, x)
    r = np.c_[x, y, axis]
    return (math.atan2(r[2, 1], r[2, 2]),
            math.atan2(-r[2, 0], math.hypot(r[2, 1], r[2, 2])),
            math.atan2(r[1, 0], r[0, 0]))


def poses(_args):
    tris = load_obj(DOME_DIR / 'base_dome.obj')
    centres = tris.mean(1)
    normals = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    areas = np.linalg.norm(normals, axis=1)
    normals = normals / np.maximum(areas, 1e-12)[:, None]

    # Centre of the dome's sphere, from the band of skin between the base ring and the top
    band = (areas > 1e-9) & (centres[:, 1] > 8) & (centres[:, 1] < 22)
    dome_centre, dome_radius = fit_sphere(centres[band])
    for _ in range(20):
        off = np.abs(np.linalg.norm(centres - dome_centre, axis=1) - dome_radius)
        dome_centre, dome_radius = fit_sphere(centres[band & (off < 0.05)])

    lines = ['<?xml version="1.0"?>',
             '<!-- Generated by tools/holo_projector_fit.py poses from base_dome.obj.',
             '     Do not edit: re-run it after re-exporting the dome. -->',
             '<robot xmlns:xacro="http://www.ros.org/wiki/xacro" >', '']
    flange = cowl_flange()
    for name, seed in SEEDS.items():
        pivot, axis, radius, gap = fit_hole(name, tris, centres, normals, areas,
                                            np.array(seed), dome_centre, flange)
        xyz = pivot * DOME_SCALE
        rpy = frame_rpy(axis)
        print(f'{name}: hole diameter {2 * radius * DOME_SCALE * 1000:.1f} mm, '
              f'flange {gap * DOME_SCALE * 1000:.1f} mm under the inner rim, '
              f'pivot {np.round(xyz, 5)}, normal {np.round(axis, 4)}')
        lines += [f'  <!-- {name}: hole diameter {2 * radius * DOME_SCALE * 1000:.1f} mm -->',
                  f'  <xacro:property name="{name}_xyz" '
                  f'value="{xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}" />',
                  f'  <xacro:property name="{name}_rpy" '
                  f'value="{rpy[0]:.6f} {rpy[1]:.6f} {rpy[2]:.6f}" />', '']
    lines += ['</robot>', '']
    (DOME_DIR / 'holo_constants.xacro').write_text('\n'.join(lines))


STL_DTYPE = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])


def read_stl(path):
    data = open(path, 'rb').read()
    count = struct.unpack('<I', data[80:84])[0]
    return np.frombuffer(data, dtype=STL_DTYPE, count=count, offset=84).copy()


def write_stl(path, mesh, offset):
    mesh['v'] = mesh['v'] - np.asarray(offset, dtype=np.float32)
    with open(path, 'wb') as f:
        f.write(b'r2_description holo projector'.ljust(80, b' '))
        f.write(struct.pack('<I', len(mesh)))
        f.write(mesh.tobytes())


def meshes(args):
    # Body: a ball with the lens barrel along +Z. Its centre is the sphere through the
    # widest part of the ball.
    body = read_stl(args.body)
    points = body['v'].reshape(-1, 3).astype(float)
    mid = (points.min(0) + points.max(0)) / 2
    r = np.linalg.norm(points[:, :2] - mid[:2], axis=1)
    ball = points[(r > 0.87 * r.max()) & (points[:, 2] < mid[2])]
    centre, radius = fit_sphere(ball)
    for _ in range(10):
        off = np.abs(np.linalg.norm(points - centre, axis=1) - radius)
        centre, radius = fit_sphere(points[off < 0.02])
    print(f'body: ball centre {np.round(centre, 3)} mm, radius {radius:.3f} mm')
    write_stl(DOME_DIR / 'holo_projector_body.stl', body, centre)

    # Cowl: a flange in the Z=0 plane with the cone rising along +Z; centre it in X/Y.
    cowl = read_stl(args.cowl)
    points = cowl['v'].reshape(-1, 3)
    offset = (points.min(0) + points.max(0)) / 2
    offset[2] = points[:, 2].min()
    print(f'cowl: flange centre {np.round(offset, 3)} mm')
    write_stl(DOME_DIR / 'holo_projector_cowl.stl', cowl, offset)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('poses').set_defaults(func=poses)
    m = sub.add_parser('meshes')
    m.add_argument('body')
    m.add_argument('cowl')
    m.set_defaults(func=meshes)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
