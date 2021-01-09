import json
from math import radians, sqrt, tan
import avlwrapper as avl

wing_span = 12
wing_aspect_ratio = 8
wing_taper = 0.3
wing_le_sweep = radians(20)
wing_dihedral = radians(4)

wing_root_chord = 2 * wing_span / (wing_aspect_ratio * (1+wing_taper))
wing_tip_chord = wing_root_chord * wing_taper

wing_root_le_pnt = avl.Point(0, 0, 0)
wing_tip_le_pnt = avl.Point(x=0.5 * wing_span * tan(wing_le_sweep),
                            y=0.5 * wing_span,
                            z=0.5 * wing_span * tan(wing_dihedral))

root_section = avl.Section(leading_edge_point=wing_root_le_pnt,
                           chord=wing_root_chord,
                           airfoil=avl.NacaAirfoil('2414'))
tip_section = avl.Section(leading_edge_point=wing_tip_le_pnt,
                          chord=wing_tip_chord,
                          airfoil=avl.NacaAirfoil('2410'))

# y_duplicate=0.0 duplicates the wing over a XZ-plane at Y=0.0
wing = avl.Surface(name='wing',
                   n_chordwise=12,
                   chord_spacing=avl.Spacing.cosine,
                   n_spanwise=20,
                   span_spacing=avl.Spacing.cosine,
                   y_duplicate=0.0,
                   sections=[root_section, tip_section])

ht_span = 5
ht_area = 6
ht_taper = 0.4
ht_sweep = radians(24)
ht_dihedral = radians(6)

ht_root_chord = 2 * ht_area / (ht_span * (1+ht_taper))
ht_tip_chord = ht_root_chord * ht_taper

elevator = avl.Control(name='elevator',
                       gain=1.0,
                       x_hinge=0.6,
                       duplicate_sign=1)

ht_root_le_pnt = avl.Point(8, 0, 0.5)
ht_tip_le_pnt = avl.Point(x=ht_root_le_pnt.x + 0.5*ht_span*tan(ht_sweep),
                          y=0.5*ht_span,
                          z=ht_root_le_pnt.z + 0.5*ht_span*tan(ht_dihedral))

root_section = avl.Section(leading_edge_point=ht_root_le_pnt,
                           chord=ht_root_chord,
                           airfoil=avl.NacaAirfoil('0012'),
                           controls=[elevator])
tip_section = avl.Section(leading_edge_point=ht_tip_le_pnt,
                          chord=ht_tip_chord,
                          airfoil=avl.NacaAirfoil('0012'),
                          controls=[elevator])
horizontal_tail = avl.Surface(name='horizontal_tail',
                              n_chordwise=12,
                              chord_spacing=avl.Spacing.cosine,
                              n_spanwise=20,
                              span_spacing=avl.Spacing.cosine,
                              y_duplicate=0.0,
                              sections=[root_section, tip_section])

mach = 0.4

wing_mac = ((2 * wing_root_chord/3) * 
            (1 + wing_taper+wing_taper**2) / 
            (1+wing_taper))

wing_area = wing_span**2 / wing_aspect_ratio

# calculate the m.a.c. leading edge location
def mac_le_pnt(root_chord, tip_chord, root_pnt, tip_pnt):
    pnt = ((2*root_chord*root_pnt[dim] +
            root_chord*tip_pnt[dim] + 
            tip_chord*root_pnt[dim] +
            2*tip_chord*tip_pnt[dim]) / 
           (3*(root_chord+tip_chord))
           for dim in range(3))
    return avl.Point(*pnt)

le_pnt = mac_le_pnt(wing_root_chord, wing_tip_chord,
                    wing_root_le_pnt, wing_tip_le_pnt)

ref_pnt = avl.Point(x=le_pnt.x + 0.25*wing_mac,
                    y=le_pnt.y, z=le_pnt.z)

acftpath = '/Users/thomasgreenhill/SUAVE/ASW22BL/ASW_22_Demo/Geometry/ASW22BL_Wing_With_Tips_Tail_Surfaces_v1.4.avl'
aircraft = avl.FileWrapper(acftpath)

 # create a session with only the geometry
session = avl.Session(geometry=aircraft)

# check if we have ghostscript
if 'gs_bin' in session.config.settings:
    img = session.save_geometry_plot()[0]
    avl.show_image(img)
else:
    session.show_geometry()

 # create a function for showing the Trefftz plot, since we'll be using it more often
def show_treffz(session):
    if 'gs_bin' in session.config.settings:
        images = session.save_trefftz_plots()
        for img in images:
            avl.show_image(img)
    else:
        for idx, _ in enumerate(session.cases):
            session.show_trefftz_plot(idx+1) # cases start from 1

simple_case = avl.Case(name='zero_aoa',
                       alpha=0)
session = avl.Session(geometry=aircraft, cases=[simple_case])

show_treffz(session)

# results are in a dictionary
result = session.run_all_cases()
print("CL = {}".format(result['zero_aoa']['Totals']['CLtot']))

session.export_run_files()

# set CL to 0.4
cl_param = avl.Parameter(name='alpha', setting='CL', value=0.4)

# trim with elevator
trim_param = avl.Parameter(name='elevator', setting='Cm', value=0.0)

trim_case = avl.Case(name='trimmed',
                     alpha=cl_param,
                     elevator=trim_param)

session = avl.Session(geometry=aircraft, cases=[trim_case])

show_treffz(session)

print("Parameters:")
print(", ".join(avl.Case.CASE_PARAMETERS.keys()))

print("\nSettings:")
print(", ".join(avl.Case.VALID_SETTINGS))

print("\nStates:")
print(", ".join(avl.Case.CASE_STATES.keys()))

base_case = avl.Case(name='sweep')

alphas = list(range(-6, 7, 2))
elevators = list(range(-15, 16, 3))
all_cases = avl.create_sweep_cases(base_case=base_case,
                                   parameters=[{'name': 'alpha',
                                                'values': alphas},
                                               {'name': 'elevator',
                                                'values': elevators}])
session = avl.Session(geometry=aircraft, cases=all_cases)

# THIS WILL FAIL!
try:
    session.run_all_cases()
except avl.session.InputError:
    print("OOPS!")

partitions = avl.partitioned_cases(all_cases)
results = {}
for partition in partitions:
    session = avl.Session(geometry=aircraft, cases=partition)
    results.update(session.run_all_cases())
    
# Write everything to json
with open('all_cases.json', 'w') as f:
    f.write(json.dumps(results))