# r2_description

The droid's URDF (xacro) and meshes. `urdf/robot_core.xacro` is the entry point; the full
robot, with its controllers, is `r2_bringup/urdf/robot.urdf.xacro`.

Joints moved by a hobby servo carry `<astromech:servo display="..."/>` (the astromech URDF
extension, see `astromech_msgs`). The servo device and the simulator find their servos by
that tag.

## Dome

`urdf/dome.xacro`: the dome, its pie panels and lower panels, and the holo projectors.
Panel hinge positions are in `meshes/dome/dome_constants.xacro`.

### Holo projectors

`urdf/holo_projector.xacro` has one macro, `holo_projector`, used three times: `hp_front`,
`hp_top` and `hp_rear`. Each is a cowl fixed to the dome and a ball-jointed body. Two servos
push and pull the ball in X and Y, so it is two revolute joints whose axes cross at the
ball's centre:

| Joint | Moves the lens | Range |
|---|---|---|
| `<name>_x_joint` | left/right | ±25° (`travel`), 0 centred |
| `<name>_y_joint` | up/down | ±25° (`travel`), 0 centred |

The projector's frame (link `<name>`) is at the ball's centre: +Z out along the lens at rest
(the dome skin's normal at the hole),
+X to the viewer's right when looking into the lens, +Y towards the top of the dome. A
positive joint position moves the lens towards + on its axis. Seen from above, +X swings
every lens anticlockwise: the front one to the droid's left, the rear one to the droid's
right. +Y tilts the front and rear ones up and the top one towards the dome's apex.

The poses are not typed in. `tools/holo_projector_fit.py poses` finds the three 65 mm holes
in `meshes/dome/base_dome.obj` and writes `meshes/dome/holo_constants.xacro`. The cowl's flange
is inside the dome, flat against the inner skin, and its cone comes out through the hole.
The holes are not cut quite square to the skin (3° to 6° off), so the lens axis is the skin's
normal, from the plane of the hole's inner rim, not the cut's axis. The pivot is under the
middle of the hole, level with the flange's base; the flange's size is read from the cowl
mesh. Re-run it after re-exporting the dome (needs numpy):

    python3 tools/holo_projector_fit.py poses

`SEEDS` at the top of the script says roughly where each hole is. A seed within 10 mm of the
hole gives the same answer; further off, or on something that is not a round hole, the
script stops and says so.

The meshes are the CSL HoloEye body and upper cowl, in millimetres.
`tools/holo_projector_fit.py meshes BODY.stl COWL.stl` re-origins the original STLs: the
body to the centre of its ball, the cowl to the centre of its flange's base, lens along +Z.

To add a projector somewhere else, call the macro with any `parent`, `xyz` and `rpy`.
